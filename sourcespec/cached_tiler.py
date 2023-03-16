"""
CachedTiler class for Catyopy.

:copyright:
    2018-2023 Claudio Satriano <satriano@ipgp.fr>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)

Adapted from
https://github.com/SciTools/cartopy/issues/732#issuecomment-191423035
"""
import os
import types
import requests
import PIL
import logging
logger = logging.getLogger(__name__.split('.')[-1])
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)


class CachedTiler(object):
    """A Cached Tiler for Cartopy."""

    def __init__(self, tiler, cache_dir=None):
        """Init with a Cartopy tiler."""
        self.tiler = tiler
        self.cache_dir = cache_dir

    def __getattr__(self, name):
        """__getattr__ method."""
        # Mimic the tiler interface, but for methods, ensure that the "self"
        # that is passed through continues to be CachedTiler, and not the
        # contained tiler instance.
        attr = getattr(self.tiler, name, None)
        if isinstance(attr, types.MethodType):
            attr = types.MethodType(attr.__func__, self)
        return attr

    def get_image(self, tile):
        """Only download a tile if it is not cached."""
        tileset_name = f'{self.tiler.__class__.__name__.lower()}'
        cache_dir = self.cache_dir
        if cache_dir is None:
            cache_dir = os.path.expanduser(
                os.path.join(
                    '~/', '.local/share/cartopy/image_tiles', tileset_name))
        else:
            cache_dir = os.path.join(cache_dir, tileset_name)
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        tile_fname = os.path.join(
            cache_dir, '_'.join(str(v) for v in tile) + '.png')
        if not os.path.exists(tile_fname):
            try:
                response = requests.get(self._image_url(tile), stream=True)
            except requests.exceptions.RequestException:
                logger.warning(
                    'Cannot download map tiles. Check internet connection.')
                # we call the original tiler, which will return an empty image
                return self.tiler.get_image(tile)
            with open(tile_fname, 'wb') as fh:
                for chunk in response:
                    fh.write(chunk)
        with open(tile_fname, 'rb') as fh:
            img = PIL.Image.open(fh)
            img = img.convert(self.desired_tile_form)
        return img, self.tileextent(tile), 'lower'
