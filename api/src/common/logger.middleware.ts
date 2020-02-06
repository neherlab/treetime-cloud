import expressWinston from 'express-winston'
import { inspect } from 'util'
import winston from 'winston'

const consoleFormat = () =>
  winston.format.printf(info => {
    const error = info?.meta?.error
    let metaString = ''
    if (error) {
      metaString = `\n\n${inspect({ ...error })}`
    }
    return `${info.timestamp} ${info.message}${metaString}`
  })

export const requestLogger = () =>
  expressWinston.logger({
    transports: [new winston.transports.Console()],
    format: winston.format.combine(
      winston.format.colorize(),
      winston.format.timestamp(),
      winston.format.align(),
      consoleFormat(),
    ),
    msg:
      'HTTP {{res.statusCode}} | {{res.responseTime}}ms | {{req.method}} | {{req.url}}',
    expressFormat: false,
    colorize: true,
  })

export const errorLogger = () =>
  expressWinston.errorLogger({
    transports: [new winston.transports.Console()],
    format: winston.format.combine(
      winston.format.colorize(),
      winston.format.timestamp(),
      winston.format.align(),
      consoleFormat(),
    ),
  })
