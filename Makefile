TILTFILE_PROD='infra/local/production.tiltfile'

.PHONY: all

all: prod

prod-up:
	tilt up -f '${TILTFILE_PROD}' --stream

prod-down:
	tilt down -f '${TILTFILE_PROD}'

prod:
	# Run `tilt up` and ensure `tilt down` runs on Ctrl+C or error
	# Probably only works on *nix platforms.
	bash -c "trap '$(MAKE) prod-down' SIGINT SIGTERM ERR; $(MAKE) prod-up"
