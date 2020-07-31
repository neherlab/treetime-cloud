FROM ubuntu:bionic-20200112

ENV DEBIAN_FRONTEND=noninteractive

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
ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"
RUN set -x \
  && curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

# Install poetry
# https://github.com/python-poetry/poetry
ENV POETRY_VERSION="1.0.3"
ENV POETRY_INSTALLER_URL="https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py"
ENV POETRY_HOME="/home/.poetry"
ENV POETRY_CACHE="/code/.cache/pypoetry"
ENV PIP_CACHE="/code/.cache/pip"
ENV PATH="${POETRY_HOME}/bin:$PATH"
RUN set -x \
  && curl -fsSL "${POETRY_INSTALLER_URL}" | /usr/bin/python3 - --version "${POETRY_VERSION}" \
  && chmod 0777 "${POETRY_HOME}/bin/poetry" \
  && sed -i "s/python/python3/g" "${POETRY_HOME}/bin/poetry"

COPY ./tools/wait-file.sh /tools/

# Install yarn
# https://github.com/yarnpkg/yarn
RUN set -x \
  && npm install -g nodemon yarn@1.22.0

ENV PYLINTHOME="/code/services/worker/.cache/pylint"

WORKDIR /code/services/worker

CMD ["poetry", "run", "python3"]
