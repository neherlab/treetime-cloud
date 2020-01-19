import '../config/dotenv'

import 'reflect-metadata'

import Server from './server/Server'

let server = new Server()
server.listen()

// This will hot-reload the Server on changes in dev mode
declare const module: NodeHotModule
if (module.hot) {
  module.hot.accept('./server/Server', async () => {
    const ServerNew = (await import('./server/Server')).default
    server.close()
    server = new ServerNew()
    server.listen()
  })
}

process.on('beforeExit', () => server.close())
