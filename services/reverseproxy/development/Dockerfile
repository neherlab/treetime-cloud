FROM nginx:1.17.8-alpine as builder

ENV UID           1000
ENV GID           1000
ENV USER          user
ENV GROUP         user

RUN set -x \
  && apk add --no-cache --virtual .build-deps \
  curl \
  inotify-tools

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

USER $USER

CMD ["nginx", "-g", "daemon off;"]
