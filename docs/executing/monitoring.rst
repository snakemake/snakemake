.. _monitoring:

========================
Monitoring
========================

`Panoptes <https://github.com/panoptes-organization/panoptes>`_ is a server (under development) that lets you monitor the execution of snakemake workflows.
Snakemake can communicate with panoptes via the ``--wms-monitor`` flag. The flag specifies the ip and port where panoptes is running (e.g. ``--wms-monitor http://127.0.0.1:5000``).

- When a snakemake pipeline is triggered with the ``--wms-monitor`` flag, an http handshake between snakemake and panoptes is performed via the panoptes ``/create_workflow`` endpoint. This returns a unique name (str(uuid.uuid4())) to snakemake for the specific pipeline. An autoincremented id is also stored on panoptes.
- All pipelines can be viewed via the ``/workflows`` endpoint or via the ``/api/workflows`` endpoint in json format.
- The status of a specific workflow can be viewed via the ``/workflow/<workflow id>`` endpoint or via the ``/api/workflow/<workflow id>`` in json format
- Similar endpoints are provided for independent jobs/rules via ``/workflow/<workflow id>/job/<job id>`` or ``/api/workflow/<workflow id>/job/<job id>``
- All jobs/rules of a specific workflow are accessible via ``/workflow/<workflow id>`` or ``/api/workflow/<workflow id>/jobs``
- The status of a workflow or a job is updated via the ``/update_workflow_status`` endpoint