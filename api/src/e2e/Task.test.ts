import 'reflect-metadata'

import path from 'path'

import request from 'supertest'

import { INestApplication } from '@nestjs/common'
import { Test } from '@nestjs/testing'

import { AppModule } from '../App.module'

import { findModuleRoot } from '../../lib/findModuleRoot'

const { moduleRoot } = findModuleRoot()

const ENDPOINT_TASK = '/api/v1/task'
const ENDPOINT_UPLOAD = '/api/v1/upload'

const DATA_PATH_BASE = path.join(moduleRoot, '..', 'examples/data/h3n2_na/h3n2_na_500') // prettier-ignore
const DATES = `${DATA_PATH_BASE}.metadata.csv`
const FASTA = `${DATA_PATH_BASE}.fasta`
const NWK = `${DATA_PATH_BASE}.nwk`

const TASK_ID = '123456'

interface FileAttachment {
  filePath: string
  fileType: string
}

interface PostParams {
  endpoint?: string
}

interface PostFileParams {
  endpoint?: string
  taskId?: string
  files: FileAttachment[]
}

describe('tasks', () => {
  let app: INestApplication | null

  beforeEach(async () => {
    const moduleRef = await Test.createTestingModule({
      imports: [AppModule],
    }).compile()

    app = moduleRef.createNestApplication()
    await app.init()
  })

  afterEach(() => {
    app?.close()
    app = null
  })

  function post<Payload>(endpoint: string, payload?: Payload) {
    const req = request(app?.getHttpServer())
      .post(endpoint)
      .set({ Accept: 'application/json' })

    if (payload) {
      req.send({ payload })
    }

    return req
  }

  function uploadFiles({ taskId, files }: PostFileParams) {
    const req = post(ENDPOINT_UPLOAD)

    req.field('taskId', taskId ?? TASK_ID)

    files.forEach(({ fileType, filePath }) => {
      req.attach(fileType, filePath)
    })

    return req
  }

  it('replies with 201/success on valid dates upload', async () => {
    const response = await uploadFiles({
      files: [{ filePath: DATES, fileType: 'DATES' }],
    })

    expect(response).toMatchObject({
      status: 201,
      body: { payload: { taskId: TASK_ID } },
    })
  })

  it('replies with 201/success on valid fasta upload', async () => {
    const response = await uploadFiles({
      files: [{ filePath: FASTA, fileType: 'FASTA' }],
    })

    expect(response).toMatchObject({
      status: 201,
      body: { payload: { taskId: TASK_ID } },
    })
  })

  it('replies with 201/success on valid nwk upload', async () => {
    const response = await uploadFiles({
      files: [{ filePath: NWK, fileType: 'NWK' }],
    })

    expect(response).toMatchObject({
      status: 201,
      body: { payload: { taskId: TASK_ID } },
    })
  })

  it('replies with 201/success on valid upload of multiple files', async () => {
    const response = await uploadFiles({
      files: [
        { filePath: FASTA, fileType: 'FASTA' },
        { filePath: DATES, fileType: 'DATES' },
        { filePath: NWK, fileType: 'NWK' },
      ],
    })

    expect(response).toMatchObject({
      status: 201,
      body: { payload: { taskId: TASK_ID } },
    })
  })

  it('replies with 201/success on valid task posted', async () => {
    const response = await post(ENDPOINT_TASK, { task: { id: TASK_ID } })

    expect(response).toMatchObject({
      status: 201,
      body: { payload: { taskId: TASK_ID } },
    })
  })
})
