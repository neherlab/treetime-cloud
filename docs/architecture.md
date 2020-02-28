## Miniservice architecture

This document briefly describes basic principles behind the miniservice
architecture and list services constituting Treetime application as a whole, as
well as tools used for development and deployment of these services.

### Miniservice architecture basics

[Miniservice architecture](https://thenewstack.io/miniservices-a-realistic-alternative-to-microservices/)
is a compromise between
[microservice architecture](https://www.martinfowler.com/articles/microservices.html)
and monolithic architecture. It allows to build modular and scalable
applications without fully investing into complexity of microservices, albeit
compromising on decoupling.

In miniservice architecture, the application consists of a set of services.
Similar to microservices, each service is a standalone application by itself,
running in a separate long-living process.

However, there are notable differences from microservices:

- The preferred means of communication in microservice world is usually a
  message queue. By contrast, miniservices use HTTP for communication where
  possible. This is simpler to get up and running, but introduces additional
  dependencies between services.

- In miniservice architecture, some related services (which would become
  different microservices in microservice architecture) are grouped together
  into one service. This reduces the complexity and management overhead, but
  increases coupling between services.

- Miniservices may share persistence and other resources. For example, a
  database is introduced as service by itself, and other services depend on it.
  By contrast, each microservice typically has its own private persistence
  capacities.

Similar to microservices, each miniservice, ideally, allows horizontal scaling.
That is, the capacity of a particular part of the application can be
independently and dynamically increased or decreased as needed, depending on
load, required performance or hosting and compute budget. This is achieved by
increasing or decreasing the number of running instances (replicas) of a
service. Scaling one service, ideally, should not affect other services.

### List of services and their functions

Our application consists of a set of services. Each service is packaged in a
container image. Containerization ensures that services can be deployed to the
same host or distributed across multiple hosts, on premises or in the cloud,
allowing to run the same application whether during local development after
installation on premises or after deployment in the cloud.

For a communication , we currently mostly use JSON APIs over HTTP transport, but
also, occasionally, AMQP protocol (message queue, for balancing tasks between
workers) and Websockets (for real-time communication between API and web
client).

Some services have multiple options for deployment, in particular, to fit
different environments: development, staging, production.

Use of containers ensures that services can be managed independently, using
simple orchestration mechanisms, such as
[docker-compose](https://docs.docker.com/compose/), as well as cloud autoscalers
and cluster offerings, for example if deployed on a
[Kubernetes](https://kubernetes.io/) cluster.

Some services that use only third party software can be replaced by cloud
offerings (e.g. AWS AmazonMQ instead of RabbitMQ, AWS load balancer instead of
Nginx, AWS S3 instead of Minio).

Below is the list of current services, the description of their functionality
and of possibilities for their deployment.

1. #### 🌍 👥 web

- frontend client application, user-facing
- allows to upload files
- allows to start a task
- shows task progress
- renders results
- allows to download results
- 🚧 current implementation: single-page application, Typescript, React, Redux,
  Saga, hosted statically with Nginx
- 🌩️deployment options:
  - served with Nginx container, locally, for development
  - served with Nginx container, on premises
  - served with Nginx container, on cloud container service (AWS ECS/EKS, Linode
    LKS)
  - any cloud object storage, e.g. AWS S3

2. #### ⚙️ 👥 api

   - RESTful API, handles requests from clients, dispatches them to appropriate
     services
   - provides unique task IDs to clients
   - accepts file uploads
   - passes file upload streams to filestore
   - accepts task requests
   - pushes tasks into the task queue
   - accepts tasks progress status from workers
   - passes task progress status to clients, using WebSockets
   - provides downloads to clients
   - ?? provides signup, login, retaining and sharing of results
   - 🚧 current implementation: Node.js, Typescript, NestJS
   - 🌩️deployment options:
     - in a container locally, for development
     - in a container, on premises
     - in a container in cloud, on AWS ECS/EKS
     - on AWS Elastic Beanstalk (cloud)

3. #### 🚀 👥 reverseproxy

   - routes requests between private and public services
   - provides SSL termination
   - provides load balancing, if necessary
   - 🚧 current implementation: Nginx
   - 🌩️deployment options:
     - Nginx container, locally, for development
     - Nginx container, on premises
     - Nginx container, on cloud container service (AWS ECS/EKS, Linode LKS)
     - any cloud reverse proxy and load balancer, e.g. AWS API Gateway, AWS Load
       Balancer

4. #### 💾 🔒 filestore

   - generic object storage
   - stores input and output files
   - accepts file uploads from API and Workers
   - 🚧 current implementation: Minio (local, on-prem) or AWS S3 (cloud)
   - 🌩️deployment options:
     - Minio container, locally, for development
     - Minio container, on premises
     - Minio container, on cloud container service (e.g. AWS ECS/EKS, Linode
       LKS)
     - any cloud object storage, e.g. AWS S3

5. #### 🤹 🔒 taskqueue

   - message queue broker to serialize task load across workers
   - API pushes tasks
   - workers pop tasks
   - 🚧 current implementation: RabbitMQ (local, on prem) or AmazonMQ (cloud)
   - 🌩️deployment:
     - RabbitMQ: in a container, locally, for development
     - RabbitMQ: in a container, on premises
     - RabbitMQ: in a container, in cloud, on AWS ECS/EKS
     - AmazonMQ, in cloud

6. #### 👷 🔒 worker

   - pops a task
   - downloads input files from filestore
   - runs computational task
   - sends progress updates to API
   - uploads output files to filestore
   - autoscaled depending on load: n persistent replicas (to reduce latency) +
     up to m dynamically allocated/deallocated replicas (to increase throughput)
   - 🚧 current implementation: Flask
   - 🌩️deployment options:
     - development: Docker - one or few persistent replicas + custom autoscaling
     - on-premises Docker + custom autoscaling
     - on-premises Kubernetes + pod autoscaling
     - cloud: AWS Batch ?? (high latency)
     - cloud: AWS ECS + AWS Fargate
     - cloud: AWS EKS + pod autoscaling

7. #### 💿 🔒 ?? database

- persist users and session info, possibly tasks
- persist user groups for sharing?
- to be used by API
- 🚧 current implementation: ?? PostgreSQL or MongoDB (local, on-prem) or Amazon
  RDS or Amazon DocumentDB (cloud)

#### Legend:

- 👥 - public service
- 🔒 - private service

### Tools and infrastructure

⚓ infra:

- infra/aws - provision AWS resources with Cloud Formation/Terraform/kubectl
- infra/compose - docker-compose files for local development and on-prem testing
- infra/k8s - manage Kubernetes cluster, if any
- infra/linode - provision Linode LKS with Ansible/Terraform/kubectl
