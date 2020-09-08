import { Controller, Get } from '@nestjs/common'

export type HealthCheckStatus = 'ok' | 'error' | 'warn'

export interface GetHealthCheckAliveResponse {
  status: HealthCheckStatus
}

export interface GetHealthCheckReadyResponse {
  status: HealthCheckStatus
}

@Controller()
export class HealthCheckController {
  @Get('/api/v1/healthcheck/alive')
  public async getHealthCheckAlive(): Promise<GetHealthCheckAliveResponse> {
    return { status: 'ok' }
  }

  @Get('/api/v1/healthcheck/ready')
  public async getHealthCheckReady(): Promise<GetHealthCheckReadyResponse> {
    return { status: 'ok' }
  }
}
