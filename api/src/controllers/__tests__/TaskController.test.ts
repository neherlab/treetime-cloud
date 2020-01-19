import 'reflect-metadata'

import path from 'path'

import express from 'express'
import request from 'supertest'

import Server from '../../server/Server'

import { findModuleRoot } from '../../../lib/findModuleRoot'

const { moduleRoot } = findModuleRoot()

let server: Server | null

let app: express.Application | null

const ENDPOINT_TASK = '/api/v1/task'
const ENDPOINT_DATES = '/api/v1/upload/dates'
const ENDPOINT_FASTA = '/api/v1/upload/fasta'
const ENDPOINT_NWK = '/api/v1/upload/nwk'

const DATA_PATH_BASE = path.join(moduleRoot, '..', 'examples/data/h3n2_na/h3n2_na_500') // prettier-ignore
const DATES = `${DATA_PATH_BASE}.metadata.csv`
const FASTA = `${DATA_PATH_BASE}.fasta`
const NWK = `${DATA_PATH_BASE}.nwk`

const TASK_ID = '123456'

beforeEach(() => {
  server = new Server({
    useRequestLogger: false,
    useErrorLogger: false,
  })

  app = server.get()
})

afterEach(() => {
  server?.close()
  server = null
  app = null
})

interface PostFileParams {
  endpoint: string
  filePath: string
  taskId?: string
}

function postFile({ endpoint, filePath, taskId }: PostFileParams) {
  return request(app)
    .post(endpoint)
    .set({ Accept: 'application/json' })
    .field('taskId', taskId ?? TASK_ID)
    .attach('file', filePath)
}

test('replies with 200/success on valid dates upload', async () => {
  const response = await postFile({
    endpoint: ENDPOINT_DATES,
    filePath: DATES,
  })

  expect(response).toMatchObject({
    status: 200,
    body: { payload: { taskId: TASK_ID } },
  })
})

test('replies with 200/success on valid fasta upload', async () => {
  const response = await postFile({
    endpoint: ENDPOINT_FASTA,
    filePath: FASTA,
  })
  expect(response).toMatchObject({
    status: 200,
    body: { payload: { taskId: TASK_ID } },
  })
})

test('replies with 200/success on valid nwk upload', async () => {
  const response = await postFile({
    endpoint: ENDPOINT_NWK,
    filePath: NWK,
    taskId: TASK_ID,
  })
  expect(response).toMatchObject({
    status: 200,
    body: { payload: { taskId: TASK_ID } },
  })
})

test('replies with 200/success on valid task posted', async () => {
  const response = await request(app)
    .post(ENDPOINT_TASK)
    .set({ Accept: 'application/json' })
    .send({ payload: { task: { id: TASK_ID } } })

  expect(response).toMatchObject({
    status: 200,
    body: { payload: { taskId: TASK_ID } },
  })
})
