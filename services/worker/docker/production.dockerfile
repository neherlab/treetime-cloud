FROM ubuntu:bionic-20200112

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV PIP_NO_CACHE_DIR=off
ENV PIP_DISABLE_PIP_VERSION_CHECK=on
ENV PIP_DEFAULT_TIMEOUT=100

RUN set -x \
  && apt-get update -qq > /dev/null \
  && apt-get install -qq --yes --force-yes --no-install-recommends \
  bash \
  build-essential \
  ca-certificates \
  coreutils \
  curl \
  fasttree \
  make \
  python3 \
  python3-pip \
  python3-venv \
  > /dev/null

# Install dockerize
# https://github.com/jwilder/dockerize
ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"
RUN set -x \
  && curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

COPY services/worker/docker/.config/pypoetry/config.toml /.config/pypoetry/

# Install poetry
# https://github.com/python-poetry/poetry
ENV POETRY_VERSION="1.0.3"
ENV POETRY_INSTALLER_URL="https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py"
ENV POETRY_VIRTUALENVS_IN_PROJECT=1
ENV POETRY_NO_INTERACTION=1
ENV POETRY_HOME="/home/.poetry"
ENV POETRY_CACHE="/code/.cache/pypoetry"
ENV VENV_PATH="/code/.cache/.venv"
ENV PIP_CACHE="/code/.cache/pip"
ENV WORKER_CACHE_DIR="/code/services/worker/.cache"
ENV WORKER_DATA_DIR="/code/services/worker/.data"
ENV PATH="${POETRY_HOME}/bin:$PATH"
RUN set -x \
  && curl -fsSL "${POETRY_INSTALLER_URL}" | /usr/bin/python3 - --version "${POETRY_VERSION}" \
  && chmod 0777 "${POETRY_HOME}/bin/poetry"

RUN set -x \
  && mkdir -p /.local /.cache $POETRY_HOME $POETRY_CACHE $PIP_CACHE $WORKER_CACHE_DIR $WORKER_DATA_DIR \
  && chmod 0777 -R /.local /.cache $POETRY_HOME $POETRY_CACHE $PIP_CACHE $WORKER_CACHE_DIR $WORKER_DATA_DIR

COPY services/worker /code/services/worker

WORKDIR /code/services/worker

RUN set -x \
  && poetry install --no-root --no-interaction

CMD bash -c "set -x \
  && export POETRY_HOME=${POETRY_HOME} \
  && export POETRY_VIRTUALENVS_IN_PROJECT=${POETRY_VIRTUALENVS_IN_PROJECT} \
  && export POETRY_NO_INTERACTION=${POETRY_NO_INTERACTION} \
  && export POETRY_CACHE=${POETRY_CACHE} \
  && export PIP_CACHE=${PIP_CACHE} \
  && export VENV_PATH=${VENV_PATH} \
  && dockerize -wait tcp://treetime-prod-taskqueue:5672 -timeout 60s >& /dev/null \
  && poetry run python src/main.py \
  "
