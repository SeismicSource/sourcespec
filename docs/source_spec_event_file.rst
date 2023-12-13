.. _source_spec_event_file:

#####################
SourceSpec Event File
#####################

The SourceSpec Event File format is a custom format, based on `YAML`_, which
allows specifying the source properties for one ore more given events.
The following properties can be specified:

*   ``event_id`` (**mandatory**): the event ID.

    .. note::
        This field must be preceded by a dash (``-``).

*   ``name`` (*optional*): the event name

*   ``hypocenter`` (**mandatory**): hypocentral location and origin time

*   ``magnitude`` (*optional*): the event magnitude (used when
    ``Mw_0_from_event_file`` in the
    :ref:`configuration_file:Configuration File` is set to ``True``)

    .. note::
        If a scalar moment or a moment tensor is given, a Mw magnitude will be
        (re)computed from it.

*   ``scalar moment`` (*optional*): event scalar moment (used when
    ``Mw_0_from_event_file`` in the
    :ref:`configuration_file:Configuration File` is set to ``True``)

    .. note::
        If a moment tensor is given, the scalar moment will be (re)computed
        from it.

*   ``focal_mechanism`` (*optional*): event focal mechanism (used for
    computing radiation pattern when the option ``rp_from_focal_mechanism`` in
    the :ref:`configuration_file:Configuration File` is set to ``True``).

    .. note::
        If a moment tensor is given, the focal mechanism will be (re)computed
        from it.

*   ``moment_tensor`` (*optional*): event moment tensor (used for computing
    focal planes and radiation pattern when the option
    ``rp_from_focal_mechanism`` in the
    :ref:`configuration_file:Configuration File` is set to ``True``;
    also used for computing moment magnitude, to be used when
    ``Mw_0_from_event_file`` is set to ``True``)

See the sample SourceSpec Event File below for more details on the
format.

Sample SourceSpec Event File
=============================

A sample SourceSpec Event File can be obtained through the command:

.. code-block:: bash

    source_spec --samplesspevent

or

.. code-block:: bash

    source_spec -y

The content of the sample SourceSpec Event File is shown below:

.. literalinclude :: ../sourcespec/config_files/ssp_event.yaml
    :language: yaml

.. _YAML: http://yaml.org/