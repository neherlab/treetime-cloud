export declare interface NodeHotLoaderWebpackPluginOptions {
  force?: boolean
  fork?: boolean | string[]
  args?: string[]
  inMemory?: boolean
  logLevel?: string
}

declare class NodeHotLoaderWebpackPlugin {
  constructor(options: NodeHotLoaderWebpackPluginOptions)
}

export default NodeHotLoaderWebpackPlugin
