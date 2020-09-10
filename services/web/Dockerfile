FROM node:12.15.0-alpine3.11 as builder

COPY services/web /code/services/web/
COPY .env.example /code/

RUN set -x \
&& cd /code \
&& cp .env.example .env

RUN set -x \
&& cd /code/services/web \
&& cp .env.experimental .env

RUN set -x \
&& cd /code/services/web \
&& yarn install --frozen-lockfile --silent

RUN set -x \
&& cd /code/services/web \
&& yarn prod:build


FROM nginx:1.17.1-alpine as runner

COPY --from=builder /code/services/web/.build/production/web /usr/share/nginx/html

# COPY nginx/nginx.conf /etc/nginx/conf.d/default.conf
