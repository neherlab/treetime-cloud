import { Request, Response } from 'express'
import expressWinston from 'express-winston'
import { inspect } from 'util'
import winston from 'winston'

interface TimedResponse extends Response {
  responseTime?: string | number
}

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
    msg(req: Request, res: Response) {
      const { url, method } = req
      const { responseTime, statusCode } = res as TimedResponse
      const timeStr = responseTime?.toString().padStart(5, ' ')
      const methodStr = method.toString().padStart(7, ' ')
      return `${statusCode} | ${timeStr}ms | ${methodStr} | ${url}`
    },
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
