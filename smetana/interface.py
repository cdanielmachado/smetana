#!/usr/bin/env python

import os
import glob
import pandas as pd
from collections import OrderedDict
from reframed import Environment
from .smetana import mip_score, mro_score, sc_score, mp_score, mu_score, minimal_environment, random_environment
from random import sample
from reframed.io.cache import ModelCache
from smetana.legacy import Community


def extract_id_from_filepath(filepath):
    filename = os.path.basename(filepath)

    if filename.endswith('.xml'):
        organism_id = filename[:-4]
    elif filename.endswith('.xml.gz'):
        organism_id = filename[:-7]
    else:
        raise IOError('Unrecognized extension in file {}. Valid extensions are .xml and .xml.gz'.format(filename))

    return organism_id


def build_cache(models, flavor=None):
    ids = [extract_id_from_filepath(model) for model in models]

    if not flavor:
        flavor = 'fbc2'

    load_args = {'flavor': flavor}

    def post_process(model):
        if 'R_ATPM' in model.reactions:
            model.reactions.R_ATPM.lb = 0

    return ModelCache(ids, models, load_args=load_args, post_processing=post_process)


def load_communities(models, communities, other, flavor):
    if len(models) == 1 and '*' in models[0]:
        pattern = models[0]
        models = glob.glob(pattern)
        if len(models) == 0:
            raise RuntimeError("No files found: {}".format(pattern))

    if other is not None:
        df = pd.read_csv(other, header=None)
        other_models = set(df[0])
    else:
        other_models = set()

    model_cache = build_cache(models, flavor)

    if communities is not None:
        df = pd.read_csv(communities, sep='\t', header=None, dtype=str)
        comm_dict = OrderedDict((name, group[1].tolist()) for name, group in df.groupby(0))
    else:
        comm_dict = {'all': model_cache.get_ids()}

    if other:
        missing = other_models - set(model_cache.get_ids())
        if missing:
            print("Warning: {} models in perturbation file missing.".format(len(missing)))

    return model_cache, comm_dict, other_models


def load_media_db(filename, sep='\t', medium_col='medium', compound_col='compound'):
    """ Load media library file. """

    data = pd.read_csv(filename, sep=sep)
    media_db = data[[medium_col, compound_col]].groupby(medium_col).agg(lambda x: list(x))

    return media_db[compound_col].to_dict()


def load_media(media, mediadb, exclude, other):
    media_db = None

    if media:
        if mediadb is None:
            raise IOError('Please provide media library with --mediadb option.')
        else:
            media_db = load_media_db(mediadb)
        media = media.split(',')
    else:
        media = [None]

    if exclude:
        df = pd.read_csv(exclude, header=None)
        excluded_mets = set(df[0])
    else:
        excluded_mets = set()

    if other:
        df = pd.read_csv(other, header=None)
        other_mets = set(df[0])
    else:
        other_mets = set()

    return media, media_db, excluded_mets, other_mets


def define_environment(medium, media_db, community, mode, verbose, min_mol_weight, use_lp):
    max_uptake = 10.0 * len(community.organisms)

    if medium:
        fmt_func = lambda x: "R_EX_M_{}_e_pool".format(x)
        env = Environment.from_compounds(media_db[medium], fmt_func=fmt_func, max_uptake=max_uptake)
        medium_id = medium
    elif mode == "global":
        env = Environment.complete(community.merged, max_uptake=max_uptake)
        medium_id = 'complete'
    elif mode == "abiotic2":
        env = random_environment(community, verbose=verbose, max_uptake=max_uptake)
        medium_id = 'random'
    else:
        env = minimal_environment(community, verbose=verbose, min_mol_weight=min_mol_weight, use_lp=use_lp,
                                  max_uptake=max_uptake)
        medium_id = "minimal"

    return medium_id, env


def run_global(comm_id, community, organisms, medium_id, excluded_mets, env, verbose, min_mol_weight, use_lp, debug):
    global_data = []
    debug_data = []

    if verbose:
        print('Running MIP for community {} on medium {}...'.format(comm_id, medium_id))

    mip, extras = mip_score(community, environment=env, verbose=verbose,
                            min_mol_weight=min_mol_weight, use_lp=use_lp, exclude=excluded_mets)

    if mip is None:
        mip = 'n/a'

    if debug and extras is not None:
        mip_ni = ','.join(sorted(extras['noninteracting_medium']))
        mip_i = ','.join(sorted(extras['interacting_medium']))
        debug_data.append((comm_id, medium_id, 'mip', 'ni', mip_ni))
        debug_data.append((comm_id, medium_id, 'mip', 'i', mip_i))

    if verbose:
        print('Running MRO for community {} on medium {}...'.format(comm_id, medium_id))

    mro, extras = mro_score(community, environment=env, verbose=verbose,
                            min_mol_weight=min_mol_weight, use_lp=use_lp, exclude=excluded_mets)

    if mro is None:
        mro = 'n/a'

    if debug and extras is not None:
        comm_medium = ','.join(sorted(extras['community_medium']))
        debug_data.append((comm_id, medium_id, 'mro', 'community', comm_medium))
        for org, values in extras['individual_media'].items():
            org_medium = ','.join(sorted(values))
            debug_data.append((comm_id, medium_id, 'mro', org, org_medium))

    global_data.append((comm_id, medium_id, len(organisms), mip, mro))

    return global_data, debug_data


def run_detailed(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight, ignore_coupling):
    smt_data = []

    exclude_bigg = {'M_{}_e'.format(x) for x in excluded_mets}

    if not ignore_coupling:
        if verbose:
            print('Running SCS for community {} on medium {}...'.format(comm_id, medium_id))

        scs = sc_score(community, environment=env, verbose=verbose)

    if verbose:
        print('Running MUS for community {} on medium {}...'.format(comm_id, medium_id))

    mus = mu_score(community, environment=env, verbose=verbose,
                   min_mol_weight=min_mol_weight)

    if verbose:
        print('Running MPS for community {} on medium {}...'.format(comm_id, medium_id))

    mps = mp_score(community, environment=env)

    pairs = [(org1, org2) for org1 in community.organisms
             for org2 in community.organisms if org1 != org2]

    for org1, org2 in pairs:
        if not ignore_coupling and scs[org1] is None:
            continue
        if mus[org1] is None:
            continue
        if mps[org2] is None:
            continue

        metabolites = (set(mus[org1]) | set(mps[org2])) - exclude_bigg

        for met in sorted(metabolites):
            mus_o1_met = mus[org1].get(met, 0)
            mps_o2_met = mps[org2].get(met, 0)

            if ignore_coupling:
                scs_o1_o2 = 'n/a'
                smt = mus_o1_met * mps_o2_met
            else:
                scs_o1_o2 = scs[org1][org2]
                smt = scs_o1_o2 * mus_o1_met * mps_o2_met
            smt_data.append((comm_id, medium_id, org1, org2, met, scs_o1_o2, mus_o1_met, mps_o2_met, smt))

    return smt_data


def run_abiotic(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight, other_mets, n, p,
                ignore_coupling):
    medium = set(env.get_compounds(fmt_func=lambda x: x[7:-7]))
    inserted = sorted(other_mets - (medium | excluded_mets))

    if len(inserted) < p:
        raise RuntimeError("Insufficient compounds ({}) to perform ({}) perturbations.".format(len(inserted), p))

    max_uptake = 10.0 * len(community.organisms)

    if n == 0:
        do_all = True
        n = len(inserted)
        if verbose:
            print('Running {} systematic abiotic perturbations with 1 compound...'.format(n))

    else:
        do_all = False
        if verbose:
            print('Running {} random abiotic perturbations with {} compounds...'.format(n, p))

    data = run_detailed(comm_id, community, medium_id, excluded_mets, env, False, min_mol_weight, ignore_coupling)

    for i in range(n):
        if do_all:
            new_compounds = list(medium) + [inserted[i]]
            new_id = "{}_{}".format(medium_id, inserted[i])
        else:
            new_compounds = list(medium) + sample(inserted, p)
            new_id = "{}_{}".format(medium_id, i + 1)

        new_env = Environment.from_compounds(new_compounds, fmt_func=lambda x: f"R_EX_M_{x}_e_pool",
                                             max_uptake=max_uptake)
        entries = run_detailed(comm_id, community, new_id, excluded_mets, new_env, False, min_mol_weight,
                               ignore_coupling)
        data.extend(entries)

    return data


def run_abiotic2(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight, n, p, ignore_coupling):
    medium = set(env.get_compounds(fmt_func=lambda x: x[7:-7]))
    removed = list(medium - excluded_mets)

    if len(removed) < p:
        raise RuntimeError("Insufficient compounds ({}) to perform ({}) perturbations.".format(len(removed), p))

    max_uptake = 10.0 * len(community.organisms)

    if n == 0:
        do_all = True
        n = len(removed)
        if verbose:
            print('Running {} systematic abiotic perturbations with 1 compound...'.format(n))
    else:
        do_all = False
        if verbose:
            print('Running {} random abiotic perturbations with {} compounds...'.format(n, p))

    data = run_detailed(comm_id, community, medium_id, excluded_mets, env, False, min_mol_weight, ignore_coupling)

    for i in range(n):
        if do_all:
            new_compounds = medium - {removed[i]}
            new_id = "{}_{}".format(medium_id, removed[i])
        else:
            new_compounds = medium - set(sample(removed, p))
            new_id = "{}_{}".format(medium_id, i + 1)

        new_env = Environment.from_compounds(new_compounds, fmt_func=lambda x: f"R_EX_M_{x}_e_pool",
                                             max_uptake=max_uptake)
        entries = run_detailed(comm_id, community, new_id, excluded_mets, new_env, False, min_mol_weight,
                               ignore_coupling)
        data.extend(entries)

    return data


def run_biotic(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight, other_models, model_cache,
               n, p, ignore_coupling):
    inserted = sorted(other_models - set(community.organisms))

    if len(inserted) < p:
        raise RuntimeError("Insufficient species ({}) to perform ({}) perturbations.".format(len(inserted), p))

    if n == 0:
        do_all = True
        n = len(inserted)
        if verbose:
            print('Running {} systematic biotic perturbations with 1 species...'.format(n))

    else:
        do_all = False
        if verbose:
            print('Running {} random biotic perturbations with {} species...'.format(n, p))

    data = run_detailed(comm_id, community, medium_id, excluded_mets, env, False, min_mol_weight, ignore_coupling)

    for i in range(n):
        if do_all:
            new_species = list(community.organisms) + [inserted[i]]
            new_id = "{}_{}".format(comm_id, inserted[i])
        else:
            new_species = list(community.organisms) + sample(inserted, p)
            new_id = "{}_{}".format(comm_id, i + 1)

        comm_models = [model_cache.get_model(org_id, reset_id=True) for org_id in new_species]
        new_community = Community(comm_id, comm_models, copy_models=False, create_biomass=False)
        entries = run_detailed(new_id, new_community, medium_id, excluded_mets, env, False, min_mol_weight,
                               ignore_coupling)
        data.extend(entries)

    return data


def export_results(mode, output, data, debug_data, zeros):
    prefix = output + '_' if output else ''

    if mode == "global":

        df = pd.DataFrame(data, columns=['community', 'medium', 'size', "mip", "mro"])
        df.to_csv(prefix + 'global.tsv', sep='\t', index=False)

        if len(debug_data) > 0:
            df = pd.DataFrame(debug_data, columns=['community', 'medium', 'key1', "key2", "data"])
            df.to_csv(prefix + 'debug.tsv', sep='\t', index=False)

    else:

        columns = ['community', 'medium', 'receiver', 'donor', 'compound', 'scs', 'mus', 'mps', 'smetana']
        df = pd.DataFrame(data, columns=columns)

        if not zeros:
            df = df.query('smetana > 0')

        df.to_csv(prefix + 'detailed.tsv', sep='\t', index=False)


def main(models, communities=None, mode=None, output=None, flavor=None, media=None, mediadb=None, zeros=False,
         verbose=False, min_mol_weight=False, use_lp=False, exclude=None, debug=False,
         other=None, n=1, p=1, ignore_coupling=False):

    other_models = other if mode == "biotic" else None
    model_cache, comm_dict, other_models = load_communities(models, communities, other_models, flavor)

    other_mets = other if mode == "abiotic" else None
    media, media_db, excluded_mets, other_mets = load_media(media, mediadb, exclude, other_mets)

    data = []
    debug_data = []

    for comm_id, organisms in list(comm_dict.items()):

        if verbose:
            print("Loading community: " + comm_id)

        comm_models = [model_cache.get_model(org_id, reset_id=True) for org_id in organisms]
        community = Community(comm_id, comm_models, copy_models=False)

        for medium in media:

            medium_id, env = define_environment(medium, media_db, community, mode, verbose, min_mol_weight, use_lp)

            if mode == "global":
                entries, debug_entries = run_global(comm_id, community, organisms, medium_id, excluded_mets, env,
                                                    verbose, min_mol_weight, use_lp, debug)

            if mode == "detailed":
                entries = run_detailed(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight,
                                       ignore_coupling)

            if mode == "abiotic":
                entries = run_abiotic(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight,
                                      other_mets, n, p, ignore_coupling)

            if mode == "abiotic2":
                entries = run_abiotic2(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight,
                                       n, p, ignore_coupling)

            if mode == "biotic":
                entries = run_biotic(comm_id, community, medium_id, excluded_mets, env, verbose, min_mol_weight,
                                     other_models, model_cache, n, p, ignore_coupling)

            data.extend(entries)

            if debug:
                debug_data.extend(debug_entries)

    export_results(mode, output, data, debug_data, zeros)

    if verbose:
        print('Done.')

