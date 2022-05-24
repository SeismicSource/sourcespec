# -*- coding: utf8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""
Save optimized PNG using PIL

:copyright:
    2022 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
import io
import PIL


def savepng(fig, figfile, **kwargs):
    """Save optimized PNG using PIL."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', **kwargs)
    buf.seek(0)
    img = PIL.Image.open(buf)
    img = img.convert('P', palette=PIL.Image.ADAPTIVE, colors=256)
    img.save(figfile, optimize=True)
    img.close()
