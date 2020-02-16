FROM python:3.6.9-slim-buster

ENV DEBIAN_FRONTEND=noninteractive

RUN set -x \
  && apt-get update -qq > /dev/null \
  && apt-get install -qq --yes --force-yes --no-install-recommends \
    build-essential \
    ca-certificates \
    coreutils \
    curl \
    fasttree \
    git \
    inotify-tools \
    make \
  > /dev/null

# Install dockerize
# https://github.com/jwilder/dockerize
ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"
RUN set -x \
  && curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

# Install poetry
# https://github.com/python-poetry/poetry
ENV POETRY_VERSION 1.0.3
ENV POETRY_INSTALLER_URL "https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py"
ENV POETRY_HOME="/home/.poetry"
ENV POETRY_CACHE="/.cache/pypoetry/cache"
ENV PIP_CACHE="/.cache/pip"
ENV PATH="${POETRY_HOME}/bin:$PATH"
RUN set -x \
  && curl -fsSL "${POETRY_INSTALLER_URL}" | python3 - --version "${POETRY_VERSION}" \
  && chmod 0777 "${POETRY_HOME}/bin/poetry" \
  && mkdir -p "${POETRY_CACHE}" && chmod -R 0777 "${POETRY_CACHE}" \
  && mkdir -p "${PIP_CACHE}" && chmod -R 0777 "${PIP_CACHE}" \
  && poetry config cache-dir "${POETRY_CACHE}" \
  && poetry config virtualenvs.create true \
  && poetry config virtualenvs.in-project true

COPY ./tools/wait-file.sh /tools/

ENV PYLINTHOME=/code/worker/.cache/pylint

WORKDIR /code/worker

CMD ["poetry", "run", "python3"]
