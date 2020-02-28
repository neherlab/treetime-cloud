# services/worker

This module implements Treetime worker instance.

This private service is used by Treetime API in order to fullfil the
computational task request from remote clients.

Treetime worker downloads input data from the filestore, processes it using
Treetime library and uploads output files back to filestore.
