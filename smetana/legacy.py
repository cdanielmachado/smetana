from collections import OrderedDict
from reframed.core.cbmodel import CBModel, CBReaction
from reframed.core.model import Compartment, Metabolite, ReactionType
from reframed.core.model import AttrOrderedDict
from warnings import warn
from copy import deepcopy
from math import inf


class CommunityNameMapping(object):
    def __init__(self, original_reaction=None, organism_reaction=None, original_metabolite=None,
                 organism_metabolite=None, extracellular_metabolite=None, community_exchange_reaction=None):
        """
        This class is used to represent mapping between original and merged community model metabolites and reactions

        Args:
            original_reaction (str): Name of reaction in original model
            organism_reaction (str): Name of reaction in merged community model
            original_metabolite (str): Name of metabolite in original model
            organism_metabolite (str): Name of metabolite in merged community model
            extracellular_metabolite (str): Name of "common environment" metabolite in merged community model
        """
        self.original_reaction = original_reaction
        self.organism_reaction = organism_reaction
        self.original_metabolite = original_metabolite
        self.organism_metabolite = organism_metabolite
        self.extracellular_metabolite = extracellular_metabolite
        self.community_exchange_reaction = community_exchange_reaction

    def __repr__(self):
        return "<orig_m: {}, org_m: {}, ex_m: {}, orig_r: {}, org_r: {}, exch_r: {}>".format(self.original_metabolite,
                                                                                             self.organism_metabolite,
                                                                                             self.extracellular_metabolite,
                                                                                             self.original_reaction,
                                                                                             self.organism_reaction,
                                                                                             self.community_exchange_reaction)


class Community(object):
    """
    This class implements a microbial community model.

    It serves as a container for multiple organisms, and can be used to merge multiple single-species models (CBModel)
    into a single multi-species model (CBModel) which is compatible with most types of constraint-based methods.
    """

    def __init__(self, community_id, models=None, copy_models=True, extracellular_compartment_id="C_e",
                 merge_extracellular_compartments=False, create_biomass=True, interacting=True,
                 exchanged_metabolites_blacklist=set()):
        """

        Args:
            community_id (str): community identifier
            models (list): list of models to be merged into single community
            copy_models (bool): If true copies for merged models  are created
            extracellular_compartment_id (str): Extracellular compartment id is used when merging extracellular compartments
            merge_extracellular_compartments (bool): Do not create organism specific extracellular compartment
            create_biomass (bool): create biomass reaction with biomass metabolites as reactants
            interacting (bool): If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
            exchanged_metabolites_blacklist (set): List of metabolites that can not be exchanged between species. This is done
             by separating 'pool' (uptake) compartment and 'pool' ('export') compartments for certain metabolites.
        """

        if not interacting and merge_extracellular_compartments:
            raise RuntimeError("Non-interacting models are not supported when merging extracellular compartment")

        self.id = community_id
        self._organisms = AttrOrderedDict()
        self._extracellular_compartment = extracellular_compartment_id  # TODO: maybe merge and compartment id arguments should be merged?
        self._merge_extracellular_compartments = merge_extracellular_compartments
        self._create_biomass = create_biomass
        self._merged_model = None
        self._copy_models = copy_models
        self._interacting = interacting
        self._organisms_exchange_reactions = {}
        self._organisms_biomass_reactions = {}
        self._exchanged_metabolites_blacklist = set(exchanged_metabolites_blacklist)

        if models is not None:
            for model in models:
                self.add_organism(model, copy_models)

    @property
    def copy_models(self):
        """
        If true copies for merged models  are created

        Returns: bool
        """
        return self._copy_models

    @property
    def create_biomass_reaction(self):
        """
        Create biomass reaction with biomass metabolites as reactants

        Returns: bool
        """
        return self._create_biomass

    @create_biomass_reaction.setter
    def create_biomass_reaction(self, value):
        """
        Create biomass reaction with biomass metabolites as reactants

        Args:
            value: bool
        """
        self._clear_merged_model()
        self._create_biomass = value

    @property
    def size(self):
        return float(len(self._organisms))

    @property
    def interacting(self):
        """
        If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink

        Returns: bool
        """
        return self._interacting

    @interacting.setter
    def interacting(self, value):
        """
        If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink

        Args:
            value: bool
        """
        self._clear_merged_model()
        self._interacting = value

    @property
    def organisms_exchange_reactions(self):
        """
        Returns dictionary containing list of reactions exchanging model metabolites with common environment.
        Dictionary keys are model ids. Values are dictionaries with keys containing exchange reaction ids and values
        containing various information about these reactions.

        Returns: dict
        """
        if not self._merged_model:
            self._merged_model = self.generate_merged_model()

        return self._organisms_exchange_reactions

    @property
    def organisms_reactions(self):
        """
        Returns dictionary containing list of community organisms specific reactions

        Returns: dict
        """
        if not self._merged_model:
            self._merged_model = self.generate_merged_model()

        return self._organisms_reactions

    @property
    def organisms_biomass_reactions(self):
        """
        Returns dictionary containing reaction exporting biomass to common environment. Keys are model ids, and values
        are reaction ids

        Returns: dict
        """
        if not self._merged_model:
            self._merged_model = self.generate_merged_model()

        return self._organisms_biomass_reactions

    @property
    def merge_extracellular_compartments(self):
        """
        Do not create organism specific extracellular compartment

        Returns: bool
        """
        return self._merge_extracellular_compartments

    @merge_extracellular_compartments.setter
    def merge_extracellular_compartments(self, value):
        """
        Do not create organism specific extracellular compartment

        Args:
            value: bool
        """
        self._clear_merged_model()
        self._merge_extracellular_compartments = value

    @property
    def merged(self):
        """
        Merged models containing every organism as separate compartment

        Returns: CBModel
        """
        if not self._merged_model:
            self._merged_model = self.generate_merged_model()

        return self._merged_model

    @property
    def organisms(self):
        """
        Dictionary of organism models which are part of the community. Keys are model ids and values are models
        Returns: dict
        """
        return self._organisms

    def __str__(self):
        return '\n'.join(self._organisms.keys())

    def _clear_merged_model(self):
        self._merged_model = None
        self._organisms_exchange_reactions = {}
        self._organisms_reactions = {}

    def add_organism(self, model, copy=True):
        """ Add an organism to this community.

        Args:
            model (CBModel): model of the organism
            copy (bool): create a copy of the given model (default: True)

        """
        self._clear_merged_model()

        if model.id in self._organisms:
            warn("Organism '{}' is already in this community".format(model.id))
        else:
            if copy:
                model = model.copy()

            self._organisms[model.id] = model

    def remove_organism(self, organism):
        """ Remove an organism from this community

        Args:
            organism (str): organism id

        """
        self._clear_merged_model()

        if organism not in self._organisms:
            warn('Organism {} is not in this community'.format(organism))
        else:
            del self._organisms[organism]

    def generate_merged_model(self):
        def _id_pattern(object_id, organism_id):
            return "{}_{}".format(object_id, organism_id)

        def _name_pattern(object_name, organism_name):
            return "{} ({})".format(object_name, organism_name)

        def _copy_object(obj, org_id, compartment=None):
            new_obj = deepcopy(obj)
            new_obj.id = _id_pattern(obj.id, org_id)
            new_obj.name = _name_pattern(obj.name, org_id)
            if compartment:
                new_obj.compartment = compartment

            return new_obj

        models_missing_extracelullar_compartment = [m.id for m in self._organisms.values()
                                                    if self._extracellular_compartment not in m.compartments]
        if models_missing_extracelullar_compartment:
            raise RuntimeError("Extracellular compartment '{}' missing from models: '{}'".format(
                self._extracellular_compartment, "', '".join(models_missing_extracelullar_compartment)))

        models_missing_biomass = [m.id for m in self._organisms.values() if not m.biomass_reaction]
        if models_missing_biomass:
            raise RuntimeError("Biomass reaction not found in models: {}".format("', '".join(models_missing_biomass)))

        merged_model = CBModel(self.id)

        organisms_biomass_metabolites = {}
        community_metabolite_exchange_lookup = {}

        for org_id, model in self._organisms.items():
            self._organisms_reactions[org_id] = []
            self._organisms_exchange_reactions[org_id] = {}
            self._organisms_biomass_reactions[org_id] = {}
            exchanged_metabolites = {m_id for r_id in model.get_exchange_reactions()
                                     for m_id in model.reactions[r_id].stoichiometry}
            #
            # Create additional extracellular compartment
            #
            if not self._merge_extracellular_compartments:
                pool_compartment = Compartment('pool', 'common pool')
                merged_model.add_compartment(pool_compartment)
                export_pool_compartment = Compartment('pool_blacklist', 'blacklisted metabolite pool')
                merged_model.add_compartment(export_pool_compartment)

            for c_id, comp in model.compartments.items():
                if c_id != self._extracellular_compartment or not self._merge_extracellular_compartments:
                    new_comp = _copy_object(comp, org_id)
                    merged_model.add_compartment(new_comp)
                elif c_id not in merged_model.compartments:
                    merged_model.add_compartment(deepcopy(comp))

            for m_id, met in model.metabolites.items():
                if met.compartment != self._extracellular_compartment or not self._merge_extracellular_compartments:
                    new_met = _copy_object(met, org_id, _id_pattern(met.compartment, org_id))
                    merged_model.add_metabolite(new_met)
                elif m_id not in merged_model.metabolites:
                    merged_model.add_metabolite(deepcopy(met))

                m_blacklisted = met.id in self._exchanged_metabolites_blacklist

                if met.id in exchanged_metabolites and not self._merge_extracellular_compartments:
                    #
                    # For blacklisted metabolites create a separate pool from which metabolites can not be reuptaken
                    #
                    if m_blacklisted and self._interacting:
                        pool_id = _id_pattern(m_id, "pool_blacklist")
                        if pool_id not in merged_model.metabolites:
                            new_met = _copy_object(met, "pool_blacklist", "pool_blacklist")
                            merged_model.add_metabolite(new_met)

                            exch_id = _id_pattern("R_EX_" + m_id, "pool_blacklist")
                            exch_name = _name_pattern(met.name, "pool (blacklist) exchange")
                            blk_rxn = CBReaction(exch_id, name=exch_name, reversible=False,
                                                 reaction_type=ReactionType.SINK)
                            blk_rxn.stoichiometry[pool_id] = -1.0
                            community_metabolite_exchange_lookup[new_met.id] = exch_id
                            merged_model.add_reaction(blk_rxn)

                    pool_id = _id_pattern(m_id, "pool")
                    if pool_id not in merged_model.metabolites:
                        new_met = _copy_object(met, "pool", "pool")
                        merged_model.add_metabolite(new_met)

                        exch_id = _id_pattern("R_EX_" + m_id, "pool")
                        exch_name = _name_pattern(met.name, "pool exchange")
                        new_rxn = CBReaction(exch_id, name=exch_name, reversible=True,
                                             reaction_type=ReactionType.EXCHANGE)
                        new_rxn.stoichiometry[pool_id] = -1.0
                        community_metabolite_exchange_lookup[new_met.id] = exch_id
                        merged_model.add_reaction(new_rxn)

            for r_id, rxn in model.reactions.items():

                is_exchange = rxn.reaction_type == ReactionType.EXCHANGE

                if not is_exchange or not self._merge_extracellular_compartments:
                    new_rxn = _copy_object(rxn, org_id)

                    for m_id, coeff in rxn.stoichiometry.items():
                        m_blacklisted = m_id in self._exchanged_metabolites_blacklist
                        if (model.metabolites[m_id].compartment != self._extracellular_compartment
                                or not self._merge_extracellular_compartments):
                            del new_rxn.stoichiometry[m_id]
                            new_id = _id_pattern(m_id, org_id)
                            new_rxn.stoichiometry[new_id] = coeff

                        if is_exchange:
                            new_rxn.reaction_type = ReactionType.OTHER
                            if (model.metabolites[m_id].compartment == self._extracellular_compartment
                                    and not self._merge_extracellular_compartments):
                                # TODO: if m_id in self._exchanged_metabolites_blacklist:
                                pool_id = _id_pattern(m_id, "pool")
                                new_rxn.stoichiometry[pool_id] = -coeff
                                cnm = CommunityNameMapping(
                                    organism_reaction=new_rxn.id,
                                    original_reaction=r_id,
                                    organism_metabolite=new_id,
                                    extracellular_metabolite=pool_id,
                                    original_metabolite=m_id,
                                    community_exchange_reaction=community_metabolite_exchange_lookup[pool_id])
                                self._organisms_exchange_reactions[org_id][new_rxn.id] = cnm

                                if not self.interacting:
                                    sink_rxn = CBReaction('Sink_{}'.format(new_id), reaction_type=ReactionType.SINK,
                                                          reversible=False)
                                    sink_rxn.stoichiometry = {new_id: -1}
                                    sink_rxn.lb = 0.0
                                    merged_model.add_reaction(sink_rxn)
                                elif m_blacklisted:
                                    pool_blacklist_id = _id_pattern(m_id, "pool_blacklist")
                                    blacklist_export_rxn = CBReaction('R_EX_BLACKLIST_{}'.format(new_id),
                                                                      reaction_type=ReactionType.OTHER,
                                                                      reversible=False)
                                    blacklist_export_rxn.stoichiometry = {new_id: -1, pool_blacklist_id: 1}
                                    blacklist_export_rxn.lb = 0.0
                                    merged_model.add_reaction(blacklist_export_rxn)

                    if is_exchange and not self._merge_extracellular_compartments:
                        new_rxn.reversible = True
                        new_rxn.lb = -inf
                        new_rxn.ub = inf if self.interacting and not m_blacklisted else 0.0

                    if rxn.id == model.biomass_reaction:
                        new_rxn.reversible = False

                    if self._create_biomass and rxn.id == model.biomass_reaction:
                        new_rxn.objective = False

                        # Add biomass metabolite to biomass equation
                        m_id = _id_pattern('Biomass', org_id)
                        name = _name_pattern('Framed biomass', org_id)
                        comp = 'pool' if not self._merge_extracellular_compartments else self._extracellular_compartment
                        biomass_met = Metabolite(m_id, name, comp)
                        merged_model.add_metabolite(biomass_met)
                        new_rxn.stoichiometry[m_id] = 1
                        organisms_biomass_metabolites[org_id] = m_id

                    self._organisms_reactions[org_id].append(new_rxn.id)
                    merged_model.add_reaction(new_rxn)

                else:
                    if is_exchange and self._merge_extracellular_compartments:
                        self._organisms_exchange_reactions[org_id][rxn.id] = CommunityNameMapping(
                            organism_reaction=r_id,
                            original_reaction=r_id,
                            extracellular_metabolite=list(rxn.stoichiometry.keys())[0],
                            original_metabolite=list(rxn.stoichiometry.keys())[0],
                            organism_metabolite=None)
                        self._organisms_reactions[org_id].append(rxn.id)

                    if r_id in merged_model.reactions:
                        continue

                    new_rxn = deepcopy(rxn)
                    new_rxn.reaction_type = ReactionType.EXCHANGE
                    if rxn.id == model.biomass_reaction and self._create_biomass:
                        new_rxn.reversible = False
                        new_rxn.objective = False

                        m_id = _id_pattern('Biomass', org_id)
                        name = _name_pattern('Biomass', org_id)
                        comp = 'pool' if not self._merge_extracellular_compartments else self._extracellular_compartment
                        biomass_met = Metabolite(m_id, name, comp)
                        merged_model.add_metabolite(biomass_met)
                        new_rxn.stoichiometry[m_id] = 1
                        organisms_biomass_metabolites[org_id] = m_id

                    merged_model.add_reaction(new_rxn)

                if r_id == model.biomass_reaction:
                    self._organisms_biomass_reactions[org_id] = new_rxn.id

        if self._create_biomass:
            biomass_rxn = CBReaction('R_Community_Growth', name="Community Growth",
                                     reversible=False, reaction_type=ReactionType.SINK, objective=1.0)
            for org_biomass in organisms_biomass_metabolites.values():
                biomass_rxn.stoichiometry[org_biomass] = -1

            merged_model.add_reaction(biomass_rxn)
            merged_model.biomass_reaction = biomass_rxn.id

        return merged_model

    def copy(self, merge_extracellular_compartments=None, copy_models=None, interacting=None, create_biomass=None,
             exchanged_metabolites_blacklist=None):
        """
        Copy model object
        Args:
            copy_models (bool): If true copies for merged models  are created
            interacting (bool): If true models will be able to exchange metabolites. Otherwise all produced metabolites will go to sink
            create_biomass (bool): create biomass reaction with biomass metabolites as reactants
        Returns:
            Community
        """
        if copy_models is None:
            copy_models = self._copy_models

        if interacting is None:
            interacting = self._interacting

        if create_biomass is None:
            create_biomass = self._create_biomass

        if merge_extracellular_compartments is None:
            merge_extracellular_compartments = self._merge_extracellular_compartments

        if exchanged_metabolites_blacklist is None:
            exchanged_metabolites_blacklist = self._exchanged_metabolites_blacklist

        copy_community = Community(self.id, models=list(self._organisms.values()),
                                   copy_models=copy_models, create_biomass=create_biomass,
                                   extracellular_compartment_id=self._extracellular_compartment,
                                   merge_extracellular_compartments=merge_extracellular_compartments,
                                   interacting=interacting,
                                   exchanged_metabolites_blacklist=exchanged_metabolites_blacklist)

        return copy_community

    def split_fluxes(self, fluxes):
        """ Decompose a flux balance solution of the merged community into organism-specific flux vectors.

        Args:
            fluxes (dict): flux distribution as a single dict

        Returns:
            dict: community flux distribution as a nested dict
        """

        comm_fluxes = OrderedDict()

        for org_id, model in self._organisms.items():
            org_fluxes = [(r_id[:-(1 + len(org_id))], val) for r_id, val in fluxes.items() if r_id.endswith(org_id)]
            comm_fluxes[org_id] = OrderedDict(org_fluxes)

        return comm_fluxes