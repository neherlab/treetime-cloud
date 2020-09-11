FROM rabbitmq:3.8.2-management


COPY "services/taskqueue/prod/etc/rabbitmq" "/etc/rabbitmq/"
COPY "services/taskqueue/prod/var/lib/rabbitmq/.erlang.cookie" "/var/lib/rabbitmq/.erlang.cookie"

RUN set -x \
  && chmod 0700 /var/lib/rabbitmq/.erlang.cookie
