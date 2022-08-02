############
Installation
############

SourceSpec requires at least Python 3.6. All the required dependencies
will be downloaded and installed during the setup process.

Using pip and PyPI (preferred method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The latest release of SourceSpec is available on the `Python Package
Index <https://pypi.org/project/sourcespec/>`__.

You can install it easily through ``pip``:

::

   pip install sourcespec

From SourceSpec GitHub releases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the latest release from the `releases
page <https://github.com/SeismicSource/sourcespec/releases>`__, in
``zip`` or ``tar.gz`` format, then:

::

   pip install sourcespec-X.Y.zip

or

::

   pip install sourcespec-X.Y.tar.gz

Where, ``X.Y`` is the version number (e.g., ``1.2``). You don’t need to
uncompress the release files yourself.

From SourceSpec GitHub repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need a recent feature that is not in the latest release (see the
``unreleased`` section in `CHANGELOG
<https://github.com/SeismicSource/sourcespec/blob/master/CHANGELOG.md>`__),
you want to use the source code from the `SourceSpec GitHub
repository <https://github.com/SeismicSource/sourcespec>`__.

For that, clone the project:

::

   git clone https://github.com/SeismicSource/sourcespec.git

(avoid using the “Download ZIP” option from the green “Code” button,
since version number is lost), then install the code from within the
``sourcespec`` main directory by running:

::

   pip install .