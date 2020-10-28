.. _monitoring:

==========
Monitoring
==========

Snakemake supports `panoptes <https://github.com/panoptes-organization/panoptes>`_ a server (under development) that lets you monitor the execution of snakemake workflows.
Snakemake communicates with panoptes via the :code:`--wms-monitor` flag. The flag specifies the ip and port where panoptes is running (e.g. :code:`--wms-monitor http://127.0.0.1:5000`).
Snakemake sends the following requests to wms monitor:

.. csv-table::
   :header: "API", "Method", "Data", "Description"
   :widths: 40, 20, 20, 60

   ":code:`/api/service-info`", "GET", "json", "Snakemake gets the status of panoptes. Snakemake continues to run if the status (:code:`json['status']`) is :code:`'running'`. In all other cases snakemake exits with an error message."
   ":code:`/create_workflow`", "GET", "json", "Snakemake gets a unique id/name :code:`str(uuid.uuid4())` for each workflow triggered."
   ":code:`/update_workflow_status`", "POST", "dictionary", "Snakemake posts updates for workflows/jobs. The dictionary sent contains the log message dictionary , the current timestamp and the unique id/name of the workflow.
   
    .. code:: python

        {
            'msg': repr(msg), 
            'timestamp': time.asctime(), 
            'id': id
        }"