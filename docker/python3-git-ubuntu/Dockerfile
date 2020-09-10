FROM ubuntu:bionic-20191029

ARG DEBIAN_FRONTEND=noninteractive

ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"

RUN apt-get update -qq \
  && apt-get install --no-install-recommends --yes --force-yes -qq \
  build-essential \
  curl \
  python3 \
  python3-dev \
  python3-pip \
  python3-setuptools \
  python3-wheel \
  > /dev/null

RUN curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

RUN set -x \
  && mkdir -p /.local /.cache /.yarn \
  && chmod 0777 /.local /.cache /.yarn

WORKDIR /code

CMD ["python3"]
