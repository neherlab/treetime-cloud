import { Module } from '@nestjs/common'

import { TaskModule } from './Task/Task.module'
import { HealthCheckModule } from './HealthCheck/HealthCheck.module'

@Module({
  imports: [TaskModule, HealthCheckModule],
})
export class AppModule {}
