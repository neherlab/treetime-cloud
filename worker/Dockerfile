FROM ubuntu:bionic-20191029

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq \
&& apt-get install --no-install-recommends --yes --force-yes -qq \
  build-essential \
  python3 \
  python3-dev \
  python3-pip \
  python3-setuptools \
  python3-wheel \
> /dev/null

COPY worker/requirements.txt /.local/

RUN mkdir -p /.local /.cache && chmod 0777 /.local /.cache

RUN pip3 install -r /.local/requirements.txt

WORKDIR "/code/worker"

CMD ["python3"]
