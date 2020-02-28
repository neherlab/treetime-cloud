# Treetime 🌲⌛

> Maximum likelihood dating and ancestral sequence inference

This is an experimental merge and remix of

- [neherlab/treetime](https://github.com/neherlab/treetime)
- [neherlab/treetime_web](https://github.com/neherlab/treetime_web)
- [neherlab/treetime_examples](https://github.com/neherlab/treetime_examples)

built with containerized miniservice architecture in mind.

## Table of contents

<!-- TOC depthFrom:2 depthTo:4 -->

- [Table of contents](#table-of-contents)
- [Getting started](#getting-started)
  - [TL;DR](#tldr)
  - [Architecture](#architecture)
  - [Directory structure](#directory-structure)
    - [Legend](#legend)
  - [Development setup](#development-setup)
    - [0. Clone git repository](#0-clone-git-repository)
    - [1. Install requirements](#1-install-requirements)
    - [2. Setup the environment](#2-setup-the-environment)
    - [3. Launch](#3-launch)
    - [4. Testing](#4-testing)
    - [5. Linting (static analysis)](#5-linting-static-analysis)
    - [6. Code style](#6-code-style)
- [Production setup](#production-setup)
- [Release setup](#release-setup)
- [Continuous integration, deployment](#continuous-integration-deployment)
- [Code of conduct](#code-of-conduct)
- [Contributing](#contributing)
- [License](#license)

<!-- /TOC -->

## Getting started

### TL;DR

If you share Kevin's "Small talk" philosophy

<a href="https://i.imgur.com/GjeQ51U.jpg">
<img src=https://i.imgur.com/GjeQ51U.jpg width="300"/>
</a>

then just copy this into your terminal window:

```bash

git clone --recursive https://github.com/ivan-aksamentov/treetime-monorepo
cd treetime-monorepo
tools/ubuntu-install-docker # if you don't have docker installed yet
cp .env{.example,}
yarn dev

```

otherwise, follow the [Development setup](#development-setup) section.

### Architecture

You can read about the architecture of the application, as well as about
included services and their functions here:
[docs/architecture](docs/architecture.md)

Take a quick look if you are a new developer!

### Directory structure

This section describes what you can find in the repository and where.

As a user you don't really need to worry about this repository, but you may find
[`examples/`](examples/) (📘) and [`docs/`](docs/) (📘) directories useful.

As a developer you are most likely interested in the actual source code (⭐) or
configuration files (🛠️).

| File or directory         | Flags  | Contents                                                                          |
| ------------------------- | ------ | --------------------------------------------------------------------------------- |
| 📁 .build/                | ⏰♻️   | Contains build artifacts                                                          |
| 📁 .cache/                | ⏰♻️   | Cache or temporary files, for example package manager cache or build system cache |
| 📁 .volumes/              | ⏰     | Writable docker volumes, to ensure persistence of files outside of containers     |
| 📁 docker/                | 🛠️⚙️   | Configuration related to docker containers and to their orchestration             |
| 📁 docs/                  | 📘     | Contains documentation                                                            |
| 📁 examples/              | 📘     | Contains examples and tutorials                                                   |
| 📁 services/api/          | ☁️⚙️⭐ | Implementation of the API service                                                 |
| 📁 services/reverseproxy/ | ☁️⚙️   | Implementation of the reverse proxy service                                       |
| 📁 services/taskqueue/    | ☁️⚙️   | Implementation of the task queue service                                          |
| 📁 services/web/          | ☁️⚙️⭐ | Implementation of the web client service                                          |
| 📁 services/worker/       | ☁️⚙️⭐ | Implementation of the worker service                                              |
| 📁 tools/                 | 🛠️     | Various developer tools and scripts                                               |
| 📄 .env                   | ⚙️     | Current environment configuration (created by the developer)                      |
| 📄 .env.example           | ⚙️     | Example environment configuration                                                 |
| 📄 package.json           | ⚙️     | Root package.json. Currently contains developer scrips                            |
| 📄 README.md              | 📘     | This document                                                                     |

#### Legend

- ⭐ - source code (start browsing here!)
- ☁️ - implements a service
- ⚙️ - contains configuration or infrastructure setup
- ⏰ - temporary
- 📘 - contains learning materials
- 📚 - implements a library
- ♻️ - safe to delete
- 🛠️ - useful tools

### Development setup

#### 0. Clone git repository

```bash
git clone --recursive https://github.com/ivan-aksamentov/treetime-monorepo
```

#### 1. Install requirements

Currently, the following are required to be installed on host machine:

- [docker](https://www.docker.com) version 19 or higher
- [docker-compose](https://docs.docker.com/compose) version 1.24 or higher
- [Node.js](https://nodejs.org/en) version 10 or higher
  ([nvm](https://github.com/nvm-sh/nvm) installer is recommended)
- [yarn](https://yarnpkg.com)

The recommended developer tools also include:

- [lazydocker](https://github.com/jesseduffield/lazydocker) - to manage
  containers from a terminal-based UI
- [ctop](https://github.com/bcicen/ctop) - to monitor performance
  characteristics of running containers

> ℹ️ As of right now, the only supported host OS is Ubuntu 18.04 LTS.

> ℹ️ Directory `tools/` contains, among others, the following utility scripts:

- [`tools/ubuntu-install-docker`](tools/ubuntu-install-docker) - installs docker
  and docker-compose, adds current user to `docker` group.

  > ⚠️ This script requires superuser privileges.

  > ⚠️ This script modifies the system. Make sure you understand what it does
  > **before** running it!

- [`tools/ubuntu-install-node`](tools/ubuntu-install-node) - installs nvm,
  Node.js and yarn
- [`tools/ubuntu-install-lazydocker`](tools/ubuntu-install-lazydocker) -
  installs lazydocker
- [`tools/ubuntu-install-ctop`](tools/ubuntu-install-ctop) - installs ctop

> ℹ️ Node.js and yarn are currently only used to run npm scripts, which further
> run services through docker-compose. Yarn scripts are likely to be replaced
> (probably with a Makefile), then the dependency on Node.js can be dropped.

#### 2. Setup the environment

Create a `.env` file with desired environment configuration. If not sure, just
copy the example one:

```bash
cp .env.example .env
```

it contains sane defaults which should be good enough to kick-start local
development. Modify as needed.

> ℹ️ Don't forget to update `.env` file every time you pull changes to the
> `.env.example`

> ⚠️ `.env` file may contain passwords, API keys and other secrets, especially
> in production environment. Keep these in a safe place!

> ⛔ Never commit your `.env` file into source control! Instead, modify
> `.env.example` file as a reference for other developers and for production
> environment.

#### 3. Launch

From project root run:

```bash
yarn dev
```

This will bring up a set of services in development mode, as described in
docker-compose file at
[`docker/docker-compose.dev.yml`](docker/docker-compose.dev.yml). In development
mode this will ensure that the necessary container images are built or pulled
from Docker Hub and that all containers are started, and their corresponding
services are built and ran.

Watch for treetime-web build, as it is the one that takes the longest. Once the
build is done, web client will be listening on `http://localhost:8080` (by
default, can be modified in `.env` file).

On code changes, most of the services will either be automatically restarted or
will have relevant modules hot reloaded. It is almost never necessary to restart
any of the services. If you find yourself in need to selectively restart a
particular service, use
[docker-compose CLI](https://docs.docker.com/compose/reference/restart/). Or,
even better, [lazydocker](https://github.com/jesseduffield/lazydocker).

Press Ctrl+C to bring all the running services down.

#### 4. Testing

TODO

#### 5. Linting (static analysis)

TODO

#### 6. Code style

TODO

## Production setup

TODO

## Release setup

TODO

## Continuous integration, deployment

TODO

## Code of conduct

TODO

## Contributing

TODO

## License

TODO
