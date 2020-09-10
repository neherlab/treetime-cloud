FROM nginx:1.17.8-alpine as builder

ENV UID           1000
ENV GID           1000
ENV USER          user
ENV GROUP         user
ENV DOCKERIZE_VERSION v0.6.1
ENV DOCKERIZE_URL "https://github.com/jwilder/dockerize/releases/download/${DOCKERIZE_VERSION}/dockerize-alpine-linux-amd64-${DOCKERIZE_VERSION}.tar.gz"

RUN set -x \
  && apk add --no-cache --virtual .build-deps \
  curl \
  inotify-tools

RUN curl -fsSL "${DOCKERIZE_URL}" | tar xfz - -C "/usr/local/bin/"

RUN set -x \
  && addgroup -g $GID -S $GROUP \
  && adduser -S -D -H -u $UID -h /home/$USER -s /sbin/nologin -G nginx -g $GROUP $USER

RUN set -x \
  && mkdir -p /home/$USER \
  && mkdir -p /web \
  && mkdir -p /sourcemaps \
  && mkdir -p /var/lib/nginx/logs \
  && mkdir -p /var/cache/nginx \
  && touch /var/lib/nginx/logs/access.log \
  && touch /var/lib/nginx/logs/error.log \
  && ln -sf /dev/stdout /var/lib/nginx/logs/access.log \
  && ln -sf /dev/stderr /var/lib/nginx/logs/error.log \
  && chown $UID:$GID -R /home/$USER \
  && chown $UID:$GID -R /var/lib/nginx \
  && chown $UID:$GID -R /var/cache/nginx \
  && chown $UID:$GID -R /web \
  && chown $UID:$GID -R /sourcemaps

COPY services/reverseproxy/production/etc/nginx /etc/nginx/

USER $USER

CMD sh -c "dockerize -wait http://treetime-prod-api:5000/api/v1/healthcheck/ready -timeout 60s >& /dev/null \
  && nginx -g 'daemon off;' \
  "
