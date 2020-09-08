#!/bin/sh

# Watches nginx config files and restarts or reloads nginx on changes.
# Useful for config development and debugging.
# 
# Requires inotify-tools:
#    apt-get update && apt-get install inotify-tools
#    apk update && apk add inotify-tools


print_colored() {
  color=$1
  message=$2
  echo -e "nginx-watch: \e[3${color}m${message}\e[0m"
}

print_error() {
  print_colored 1 "${1}"
}

print_success() {
  print_colored 2 "${1}"
}

print_info() {
  print_colored 4 "${1}"
}

# First run
nginx
sleep 1
[ -f /tmp/nginx.pid ] && read -r nginx_pid </tmp/nginx.pid
if [ -n "${nginx_pid}" ] && [ -d "/proc/${nginx_pid}" ]; then
  print_success "nginx is running with PID: $nginx_pid"
else
  print_error "nginx failed to start"
fi

# Wait for changes
while inotifywait --event modify,move,create,delete --quiet -qq --recursive '/etc/nginx'; do
    print_info "nginx configuration change detected"

    [ -f /tmp/nginx.pid ] && read -r nginx_pid </tmp/nginx.pid
    if [ -n "${nginx_pid}" ] && [ -d "/proc/${nginx_pid}" ]; then
      # If already running, test the config and try to reload
      if nginx -t >/dev/null; then
          nginx -s reload;
          print_success "nginx configuration reloaded successfully"
      else
          print_error "nginx configuration failed to reload"
          print_info "nginx running with last sucessful configuration"
      fi
    else
      # If not running, run
      echo "nginx-watch: starting nginx"
      nginx
    fi

    [ -f /tmp/nginx.pid ] && read -r nginx_pid </tmp/nginx.pid
    if [ -n "${nginx_pid}" ] && [ -d "/proc/${nginx_pid}" ]; then
      print_success "nginx is running with PID: $nginx_pid"
    else
      print_error "nginx failed to start. Last PID was: $nginx_pid"
    fi
done
