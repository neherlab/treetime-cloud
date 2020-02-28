# services/api

This module implements Treetime API server.
It is to be deployed on a cloud hosting.
 
It provides a RESTful API that:

 - isolates multiple client instances (whether web or cli)
 
 - handles file uploads and downloads
 
 - launches Treetime workers to fullfil processing tasks
 
 - provides input data to Treetime workers and provides output data back to clients
  
 - eventually, persists results
