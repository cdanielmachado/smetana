from __future__ import division
from builtins import map
from builtins import range
from past.utils import old_div
from framed.community.model import Community
from framed.experimental.medium import minimal_medium
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from framed import Environment

from collections import Counter
from itertools import combinations, chain
from warnings import warn


def species_coupling_score(community, environment=None, min_growth=0.1, n_solutions=100, verbose=True, abstol=1e-6,
                           use_pool=False):
    """
    Calculate frequency of community species dependency on each other

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): microbial community
        environment (Environment): metabolic environment (optional)
        min_growth (float): minimum growth rate (default: 0.1)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        n_solutions (int): number of alternative solutions to calculate (default: 100)

    Returns:
        dict: Keys are dependent organisms, values are dictionaries with required organism frequencies
    """

    community = community.copy(copy_models=False, interacting=True, create_biomass=False,
                               merge_extracellular_compartments=False)

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)

    for b in community.organisms_biomass_reactions.values():
        community.merged.reactions[b].lb = 0

    solver = solver_instance(community.merged)

    for org_id, rxns in community.organisms_reactions.items():
        org_var = 'y_{}'.format(org_id)
        solver.add_variable(org_var, 0, 1, vartype=VarType.BINARY, update_problem=False)

    solver.update()

    bigM = 100
    for org_id, rxns in community.organisms_reactions.items():
        org_var = 'y_{}'.format(org_id)
        for r_id in rxns:
            if r_id == community.organisms_biomass_reactions[org_id]:
                continue
            solver.add_constraint('c_{}_lb'.format(r_id), {r_id: 1, org_var: bigM}, '>', 0, update_problem=False)
            solver.add_constraint('c_{}_ub'.format(r_id), {r_id: 1, org_var: -bigM}, '<', 0, update_problem=False)

    solver.update()

    scores = {}

    for org_id, biomass_id in community.organisms_biomass_reactions.items():
        other = {o for o in community.organisms if o != org_id}
        solver.add_constraint('SMETANA_Biomass', {community.organisms_biomass_reactions[org_id]: 1}, '>', min_growth)
        objective = {"y_{}".format(o): 1.0 for o in other}


        if not use_pool:
            previous_constraints = []
            donors_list = []
            failed = False

            for i in range(n_solutions):
                sol = solver.solve(objective, minimize=True, get_values=list(objective.keys()))

                if sol.status != Status.OPTIMAL:
                    failed = i == 0
                    break

                donors = [o for o in other if sol.values["y_{}".format(o)] > abstol]
                donors_list.append(donors)

                previous_con = 'iteration_{}'.format(i)
                previous_constraints.append(previous_con)
                previous_sol = {"y_{}".format(o): 1 for o in donors}
                solver.add_constraint(previous_con, previous_sol, '<', len(previous_sol) - 1)

            solver.remove_constraints(['SMETANA_Biomass'] + previous_constraints)

            if not failed:
                donors_list_n = float(len(donors_list))
                donors_counter = Counter(chain(*donors_list))
                scores[org_id] = {o: old_div(donors_counter[o],donors_list_n) for o in other}
            else:
                if verbose:
                    warn('SCS: Failed to find a solution for growth of ' + org_id)
                scores[org_id] = None

        else:
            sols = solver.solve(objective, minimize=True, get_values=list(objective.keys()), pool_size=n_solutions, pool_gap=0.1)
            solver.remove_constraint('SMETANA_Biomass')

            if len(sols) == 0:
                scores[org_id] = None
                if verbose:
                    warn('SCS: Failed to find a solution for growth of ' + org_id)
            else:
                donor_count = [o for sol in sols for o in other if sol.values["y_{}".format(o)] > abstol]
                donor_count = Counter(donor_count)
                scores[org_id] = {o: old_div(donor_count[o], float(len(sols))) for o in other}

    return scores


def metabolite_uptake_score(community, environment=None, min_mass_weight=False, min_growth=0.1, max_uptake=10.0,
                            abstol=1e-6, validate=False, n_solutions=100, pool_gap=0.5, verbose=True,
                            exclude=None):  #TODO: implement excluded
    """
    Calculate frequency of metabolite requirement for species growth

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): microbial community
        environment (Environment): metabolic environment
        min_mass_weight (bool): Prefer smaller compounds (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        n_solutions (int): number of alternative solutions to calculate (default: 100)

    Returns:
        dict: Keys are organism names, values are dictionaries with metabolite frequencies 
        dict: Extra information
    """

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)

    scores = {}

    solver = solver_instance(community.merged)

    for org_id, exchange_rxns in community.organisms_exchange_reactions.items():
        biomass_reaction = community.organisms_biomass_reactions[org_id]
        community.merged.biomass_reaction = biomass_reaction

        medium_list, sols = minimal_medium(community.merged, exchange_reactions=list(exchange_rxns.keys()),
                                           min_mass_weight=min_mass_weight, min_growth=min_growth,
                                           n_solutions=n_solutions, max_uptake=max_uptake, validate=validate,
                                           abstol=abstol, use_pool=True, pool_gap=pool_gap, solver=solver,
                                           warnings=verbose)

        if medium_list:
            counter = Counter(chain(*medium_list))
            scores[org_id] = {cnm.original_metabolite: old_div(counter[ex], float(len(medium_list)))
                              for ex, cnm in exchange_rxns.items()}
        else:
            if verbose:
                warn('MUS: Failed to find a minimal growth medium for ' + org_id)
            scores[org_id] = None

    return scores


def metabolite_production_score(community, environment=None, abstol=1e-3, exclude=None): #TODO: implement excluded
    """
    Discover metabolites which species can produce in community

    Zelezniak A. et al, Metabolic dependencies drive species co-occurrence in diverse microbial communities (PNAS 2015)

    Args:
        community (Community): community object
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)

    Returns:
        dict: Keys are model names, values are list with produced compounds
        dict: Extra information
    """

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)
        env_compounds = environment.get_compounds(format_str="\'{}\'[5:-5]")
    else:
        env_compounds = set()

    for exchange_rxns in community.organisms_exchange_reactions.values():
        for r_id in exchange_rxns.keys():
            rxn = community.merged.reactions[r_id]
            if rxn.ub is None:
                rxn.ub = 1000

    solver = solver_instance(community.merged)

    scores = {}

    for org_id, exchange_rxns in community.organisms_exchange_reactions.items():
        scores[org_id] = {}

        remaining = [r_id for r_id, cnm in exchange_rxns.items() if cnm.original_metabolite not in env_compounds]

        while len(remaining) > 0:
            sol = solver.solve(linear={r_id: 1 for r_id in remaining}, minimize=False, get_values=remaining)

            if sol.status != Status.OPTIMAL:
                break

            blocked = [r_id for r_id in remaining if sol.values[r_id] < abstol]

            if len(blocked) == len(remaining):
                break

            for r_id in remaining:
                if sol.values[r_id] >= abstol:
                    cnm = exchange_rxns[r_id]
                    scores[org_id][cnm.original_metabolite] = 1

            remaining = blocked

        for r_id in remaining:
            cnm = exchange_rxns[r_id]
            scores[org_id][cnm.original_metabolite] = 0

    return scores


def mip_score(community, environment=None, min_mass_weight=False, min_growth=0.1, direction=-1, max_uptake=10,
              validate=False, verbose=True, exclude=None):
    """
    Implements the metabolic interaction potential (MIP) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is calculated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)

    Returns:
        float: MIP score
    """

    noninteracting = community.copy(copy_models=False, interacting=False)
    exch_reactions = set(community.merged.get_exchange_reactions())

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)
        environment.apply(noninteracting.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    interacting_medium, sol1 = minimal_medium(community.merged, direction=direction, exchange_reactions=exch_reactions,
                                              min_mass_weight=min_mass_weight, min_growth=min_growth,
                                              max_uptake=max_uptake, validate=validate, warnings=verbose)

    noninteracting_medium, sol2 = minimal_medium(noninteracting.merged, exchange_reactions=exch_reactions,
                                                 direction=direction, min_mass_weight=min_mass_weight,
                                                 min_growth=min_growth, max_uptake=max_uptake, validate=validate,
                                                 warnings=verbose)

    if noninteracting_medium is None:
        if verbose:
            warn('MIP: Failed to find a valid solution for non-interacting community')
        return None, None
    else:

        if exclude is not None:
            exclude_rxns = {'R_EX_M_{}_e_pool'.format(x) for x in exclude}
            interacting_medium = set(interacting_medium) - exclude_rxns
            noninteracting_medium = set(noninteracting_medium) - exclude_rxns

        score = len(noninteracting_medium) - len(interacting_medium)

    extras = {'noninteracting_medium': noninteracting_medium, 'interacting_medium': interacting_medium,
              'noninteracting_solution': sol2, 'interacting_solution': sol1}

    return score, extras


def mro_score(community, environment=None, direction=-1, min_mass_weight=False, min_growth=0.1, max_uptake=10,
              validate=False, verbose=True, exclude=None):
    """
    Implements the metabolic resource overlap (MRO) score as defined in (Zelezniak et al, 2015).

    Args:
        community (Community): microbial community model
        environment (Environment): Metabolic environment in which the SMETANA score is colulated
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        extracellular_id (str): extracellular compartment id
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 0.1)
        max_uptake (float): maximum uptake rate (default: 10)

    Returns:
        float: MRO score
    """

    noninteracting = community.copy(copy_models=False, interacting=False, create_biomass=True)
    exch_reactions = set(community.merged.get_exchange_reactions())

    if environment:
        environment.apply(community.merged, inplace=True, warning=False)
        environment.apply(noninteracting.merged, inplace=True, warning=False)
        exch_reactions &= set(environment)

    noninteracting_medium, sol = minimal_medium(noninteracting.merged, exchange_reactions=exch_reactions,
                                                direction=direction, min_mass_weight=min_mass_weight,
                                                min_growth=min_growth, max_uptake=max_uptake, validate=validate,
                                                warnings=verbose)

    solutions = [sol]

    if sol.status != Status.OPTIMAL:
        if verbose:
            warn('MRO: Failed to find a valid solution for non-interacting community')
        return None, None

    # anabiotic environment is limited to non-interacting community minimal media
    noninteracting_exch = set(noninteracting_medium)
    noninteracting_env = Environment.from_reactions(noninteracting_exch, max_uptake=max_uptake)
    noninteracting_env.apply(community.merged, inplace=True)

    individual_media = {}

    if exclude is not None:
        exclude = {'M_{}_e'.format(x) for x in exclude}
    else:
        exclude = {}

    solver = solver_instance(community.merged)
    for org_id in community.organisms:
        biomass_reaction = community.organisms_biomass_reactions[org_id]
        community.merged.biomass_reaction = biomass_reaction

        org_noninteracting_exch = community.organisms_exchange_reactions[org_id]

        medium, sol = minimal_medium(community.merged, exchange_reactions=org_noninteracting_exch, direction=direction,
                                     min_mass_weight=min_mass_weight, min_growth=min_growth, max_uptake=max_uptake,
                                     validate=validate, solver=solver, warnings=verbose)
        solutions.append(sol)

        if sol.status != Status.OPTIMAL:
            warn('MRO: Failed to find a valid solution for: ' + org_id)
            return None, None

        individual_media[org_id] = {org_noninteracting_exch[r].original_metabolite for r in medium} - exclude

    pairwise = {(o1, o2): individual_media[o1] & individual_media[o2] for o1, o2 in combinations(community.organisms, 2)}

    numerator = len(individual_media) * sum(map(len, pairwise.values()))
    denominator = float(len(pairwise) * sum(map(len, individual_media.values())))

    score = old_div(numerator, denominator) if denominator != 0 else None
    extras = {'noninteracting_medium': noninteracting_medium, 'individual_media': individual_media, 
              'pairwise': pairwise, 'solutions': solutions}

    return score, extras
