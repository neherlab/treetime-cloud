import { Module } from '@nestjs/common'

import { HealthCheckController } from './HealthCheck.controller'

@Module({
  controllers: [HealthCheckController],
})
export class HealthCheckModule {}
