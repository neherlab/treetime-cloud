FROM alpine:3.11 as builder

ENV NGINX_VERSION 1.17.8

RUN set -x \
&& apk add --no-cache --virtual .build-deps \
  curl \
  gcc \
  gettext \
  git \
  gnupg \
  inotify-tools \
  libc-dev \
  libtool \
  linux-headers \
  make \
  openssl-dev \
  pcre-dev \
  zlib-dev


RUN set -x \
&& mkdir "/src" \
&& cd "/src" \
&& curl -OJfsSL "https://nginx.org/download/nginx-${NGINX_VERSION}.tar.gz" \
&& curl -OJfsSL "https://nginx.org/download/nginx-${NGINX_VERSION}.tar.gz.asc" \
&& curl -OJfsSL "https://nginx.org/keys/mdounin.key" \
&& gpg --import "mdounin.key" \
&& gpg --verify "nginx-${NGINX_VERSION}.tar.gz.asc" "nginx-${NGINX_VERSION}.tar.gz"


RUN set -x \
&& cd "/src" \
&& tar -xf "nginx-${NGINX_VERSION}.tar.gz"


RUN set -x \
&& cd "/src" \
&& git clone --recursive --depth=1 "https://github.com/google/ngx_brotli" \
&& cd ngx_brotli \
&& git checkout e505dce68acc190cc5a1e780a3b0275e39f160ca


RUN set -x \
&& cd "/src/nginx-${NGINX_VERSION}" \
&& ./configure \
  --prefix=/var/lib/nginx \
  --sbin-path=/usr/sbin/nginx \
  --conf-path=/etc/nginx/nginx.conf \
  --pid-path=/run/nginx/nginx.pid \
  --lock-path=/run/nginx/nginx.lock \
  --http-client-body-temp-path=/var/lib/nginx/tmp/client_body \
  --http-proxy-temp-path=/var/lib/nginx/tmp/proxy \
  --http-fastcgi-temp-path=/var/lib/nginx/tmp/fastcgi \
  --http-uwsgi-temp-path=/var/lib/nginx/tmp/uwsgi \
  --http-scgi-temp-path=/var/lib/nginx/tmp/scgi \
  --with-perl_modules_path=/usr/lib/perl5/vendor_perl \
  --user="nginx" \
  --group="nginx" \
  --add-module="../ngx_brotli" \
  --with-file-aio \
  --with-threads \
  --with-http_gunzip_module \
  --with-http_gzip_static_module \
  --with-http_realip_module \
  --with-http_ssl_module \
  --with-http_v2_module \
  --without-http_auth_basic_module \
  --without-http_autoindex_module \
  --without-http_browser_module \
  --without-http_empty_gif_module \
  --without-http_fastcgi_module \
  --without-http_geo_module \
  --without-http_grpc_module \
  --without-http_memcached_module \
  --without-http_mirror_module \
  --without-http_scgi_module \
  --without-http_split_clients_module \
  --without-http_ssi_module \
  --without-http_userid_module \
  --without-http_uwsgi_module \
  --without-mail_imap_module \
  --without-mail_pop3_module \
  --without-mail_smtp_module \
  --without-poll_module \
  --without-select_module \
  --without-stream_access_module \
  --without-stream_geo_module \
  --without-stream_limit_conn_module \
  --without-stream_map_module \
  --without-stream_return_module \
  --without-stream_split_clients_module \
  --without-stream_upstream_hash_module \
  --without-stream_upstream_least_conn_module \
  --without-stream_upstream_zone_module


RUN set -x \
&& cd "/src/nginx-${NGINX_VERSION}" \
&& make -j$(getconf _NPROCESSORS_ONLN) \
&& make install \
&& strip /usr/sbin/nginx




FROM alpine:3.11 as runner


ENV UID           1000
ENV GID           1000


COPY --from=builder /usr/bin/envsubst /usr/bin/envsubst
COPY --from=builder /usr/bin/inotifywait /usr/bin/inotifywait
COPY --from=builder /usr/lib/libinotifytools.so.0 /usr/lib/libinotifytools.so.0
COPY --from=builder /usr/sbin/nginx /usr/sbin/nginx

RUN set -x \
&& apk add --no-cache libcrypto1.1 libssl1.1 musl zlib pcre libgcc


RUN set -x \
&& addgroup -g $GID -S nginx \
&& adduser -S -D -H -u $UID -h /var/cache/nginx -s /sbin/nologin -G nginx -g nginx nginx


RUN set -x \
&& mkdir -p /web \
&& mkdir -p /sourcemaps \
&& mkdir -p /var/lib/nginx/logs \
&& touch /var/lib/nginx/logs/access.log \
&& touch /var/lib/nginx/logs/error.log \
&& ln -sf /dev/stdout /var/lib/nginx/logs/access.log \
&& ln -sf /dev/stderr /var/lib/nginx/logs/error.log \
&& chown $UID:$GID -R /var/lib/nginx \
&& chown $UID:$GID -R /web \
&& chown $UID:$GID -R /sourcemaps

USER nginx



EXPOSE 8080
EXPOSE 44443

STOPSIGNAL SIGTERM

CMD ["nginx", "-g", "daemon off;"]
