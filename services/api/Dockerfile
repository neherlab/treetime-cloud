FROM node:12.15.0-alpine3.11 as builder

ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-alpine-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"

RUN apk add --update --no-cache \
  curl \
  > /dev/null

RUN curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

COPY services/api /code/services/api/
COPY .env.example /code/

RUN set -x \
  && cd /code \
  && cp .env.example .env \
  && cd /code/services/api \
  && yarn install --frozen-lockfile --silent \
  && yarn prod:build


FROM node:12.15.0-alpine3.11 as runner

COPY --from=builder /code/services/api/.build/production/node/api.js /code/services/api/.build/production/node/api.js
COPY --from=builder /code/services/api/package.json /code/services/api/package.json
COPY --from=builder /code/.env /code/.env
COPY --from=builder /usr/local/bin/dockerize /usr/local/bin/

RUN set -x \
  && cd /code/services/api \
  && yarn install --frozen-lockfile --silent --production

ENTRYPOINT ["docker-entrypoint.sh"]

WORKDIR /code/services/api

CMD sh -c "dockerize -wait tcp://treetime-prod-taskqueue:5672 -timeout 60s >& /dev/null \
  && node .build/production/node/api.js \
  "
