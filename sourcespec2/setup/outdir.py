# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Setup functions for sourcespec.

:copyright:
    2012 Claudio Satriano <satriano@ipgp.fr>

    2013-2014 Claudio Satriano <satriano@ipgp.fr>,
              Emanuela Matrullo <matrullo@geologie.ens.fr>,
              Agnes Chounet <chounet@ipgp.fr>

    2015-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import os
import shutil
import contextlib
from .config import config


def save_config(eventid=None):
    """Save config file to output dir."""
    if eventid is None:
        try:
            eventid = config.event.event_id
        except AttributeError as err:
            raise RuntimeError(
                'No event ID found. Cannot save config file.'
            ) from err
    # Actually, it renames the file already existing.
    src = os.path.join(config.options.outdir, 'source_spec.conf')
    dst = os.path.join(config.options.outdir, f'{eventid}.ssp.conf')
    # On Windows, dst file must not exist
    with contextlib.suppress(Exception):
        os.remove(dst)
    os.rename(src, dst)


def get_outdir_path(eventid=None):
    """Construct full path to output directory"""
    if eventid is None:
        try:
            eventid = config.event.event_id
        except AttributeError as err:
            raise RuntimeError(
                'No event ID found. Cannot create output directory.'
            ) from err
    src = config.options.outdir
    run_id = config.options.run_id
    run_id_subdir = config.options.run_id_subdir
    # TODO: does next line also work if no tmpdir has been created first?
    outdir = os.path.split(src)[0]
    outdir = os.path.join(outdir, str(eventid))
    if run_id and run_id_subdir:
        outdir = os.path.join(outdir, str(run_id))
    elif run_id:
        outdir += f'_{run_id}'

    return outdir


def move_outdir(eventid=None):
    """Move outdir to a new dir named from evid (and optional run_id)."""
    src = config.options.outdir
    dst = get_outdir_path(eventid)
    # Create destination
    if not os.path.exists(dst):
        os.makedirs(dst)
    # Copy all files into destination
    file_names = os.listdir(src)
    for file_name in file_names:
        shutil.copyfile(
            os.path.join(src, file_name),
            os.path.join(dst, file_name)
        )
    # Old outdir cannot be removed yet, because the log file is still opened
    config.options.oldoutdir = src
    config.options.outdir = dst


def remove_old_outdir():
    """Try to remove the old outdir."""
    try:
        oldoutdir = config.options.oldoutdir
        shutil.rmtree(oldoutdir)
    except Exception:
        return
