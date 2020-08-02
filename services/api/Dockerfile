FROM node:12.15.0-alpine3.11 as builder

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

RUN set -x \
&& cd /code/services/api \
&& yarn install --frozen-lockfile --silent --production

ENTRYPOINT ["docker-entrypoint.sh"]

WORKDIR /code/services/api

CMD ["node", ".build/production/node/api.js"]
