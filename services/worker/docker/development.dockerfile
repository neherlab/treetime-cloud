FROM ubuntu:bionic-20200112

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV DOCKERIZE_VERSION=v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"
ENV POETRY_VERSION="1.0.3"
ENV POETRY_INSTALLER_URL="https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py"
ENV POETRY_HOME="/home/.poetry"
ENV POETRY_CACHE="/code/.cache/pypoetry"
ENV VENV_PATH="/code/.cache/.venv"
ENV PIP_CACHE="/code/.cache/pip"
ENV PIP_NO_CACHE_DIR=off
ENV PIP_DISABLE_PIP_VERSION_CHECK=on
ENV PIP_DEFAULT_TIMEOUT=100
ENV WORKER_CACHE_DIR="/code/services/worker/.cache"
ENV WORKER_DATA_DIR="/code/services/worker/.data"
ENV PATH="${POETRY_HOME}/bin:$PATH"
ENV PYLINTHOME="/code/services/worker/.cache/pylint"

RUN set -x \
  && apt-get update -qq > /dev/null \
  && apt-get install -qq --yes --force-yes --no-install-recommends \
  bash \
  build-essential \
  ca-certificates \
  coreutils \
  curl \
  fasttree \
  git \
  inotify-tools \
  make \
  nodejs \
  npm \
  python3 \
  python3-pip \
  python3-venv \
  > /dev/null

# Install dockerize
# https://github.com/jwilder/dockerize
RUN set -x \
  && curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

# Install poetry
# https://github.com/python-poetry/poetry
RUN set -x \
  && curl -fsSL "${POETRY_INSTALLER_URL}" | /usr/bin/python3 - --version "${POETRY_VERSION}" \
  && chmod 0777 "${POETRY_HOME}/bin/poetry" \
  && sed -i "s/python/python3/g" "${POETRY_HOME}/bin/poetry"

# Create required directories and give full permissions
# TODO: create a nonroot user, chown and narrow down the permissions instead
RUN set -x \
  && mkdir -p /.local /.cache /.yarn $POETRY_HOME $POETRY_CACHE $PIP_CACHE $WORKER_CACHE_DIR $WORKER_DATA_DIR \
  && chmod 0777 -R /.local /.cache /.yarn $POETRY_HOME $POETRY_CACHE $PIP_CACHE $WORKER_CACHE_DIR $WORKER_DATA_DIR

# Install yarn
# https://github.com/yarnpkg/yarn
RUN set -x \
  && npm install -g nodemon yarn@1.22.0

WORKDIR /code/services/worker

CMD ["poetry", "run", "python3"]
