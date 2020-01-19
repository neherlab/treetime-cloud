import expressWinston from 'express-winston'
import winston from 'winston'

const consoleFormat = () =>
  winston.format.printf(info => `${info.timestamp} ${info.message}`)

export const requestLogger = () =>
  expressWinston.logger({
    transports: [new winston.transports.Console()],
    format: winston.format.combine(
      winston.format.colorize(),
      winston.format.timestamp({ format: 'YYYY.MM.DD HH:mm:ss.SSS' }),
      winston.format.align(),
      consoleFormat(),
    ),
    meta: false,
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
      winston.format.timestamp({ format: 'YYYY.MM.DD HH:mm:ss.SSS' }),
      winston.format.align(),
      consoleFormat(),
    ),
  })
