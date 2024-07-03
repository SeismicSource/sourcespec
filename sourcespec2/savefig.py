# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Save Matplotlib figure. Optimize PNG format using PIL.

:copyright:
    2022-2024 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import io
import logging
import warnings
from PIL import Image
# Reduce logging level for PIL to avoid DEBUG messages
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.WARNING)
# Reduce logging level for fontTools to avoid DEBUG messages
mpl_logger = logging.getLogger('fontTools')
mpl_logger.setLevel(logging.WARNING)
# Silence PIL warnings about transparency needing RGBA, since we remove the
# alpha channel anyway (see below)
warnings.filterwarnings('ignore', message='Palette images with Transparency')


def savefig(fig, figfile, fmt, quantize_colors=True, **kwargs):
    """Save Matplotlib figure. Optimize PNG format using PIL."""
    if fmt == 'png':
        buf = io.BytesIO()
        fig.savefig(buf, format='png', **kwargs)
        buf.seek(0)
        img = Image.open(buf)
        if quantize_colors:
            # pylint: disable=maybe-no-member
            img = img.convert('P', palette=Image.ADAPTIVE, colors=256)
        else:
            # just remove the alpha channel
            img = img.convert('RGB')
        img.save(figfile, optimize=True)
        img.close()
    else:
        fig.savefig(figfile, **kwargs)
