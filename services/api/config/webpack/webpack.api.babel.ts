import '../dotenv'

import path from 'path'

import ExtraWatchWebpackPlugin from 'extra-watch-webpack-plugin'
import webpack from 'webpack'
import nodeExternals from 'webpack-node-externals'
import StartServerPlugin from 'start-server-webpack-plugin'

import webpackLoadJavascript from './lib/webpackLoadJavascript'
import webpackTerser from './lib/webpackTerser'
import webpackTsChecker from './lib/webpackTsChecker'

import { findModuleRoot } from '../../lib/findModuleRoot'
import { getenv } from '../../lib/getenv'

import babelConfig from '../../babel.config'

const MODE: 'development' | 'production' = getenv('NODE_ENV') === 'development' ? 'development' : 'production' // prettier-ignore

const production = MODE === 'production'
const development = MODE === 'development'
const sourceMaps = true

const { moduleRoot, pkg } = findModuleRoot()
const buildPath = path.join(moduleRoot, '.build', MODE)
const outputFilename = 'api.js' as const

function alias(development: boolean) {
  let productionAliases = {}

  if (development) {
    productionAliases = {
      ...productionAliases,
    }
  }

  return productionAliases
}

module.exports = {
  mode: MODE,
  bail: true,
  name: 'api',
  target: 'node',
  devtool: 'cheap-module-source-map',
  stats: {
    all: false,
    errors: true,
    warnings: true,
    moduleTrace: true,
    colors: true,
  },
  performance: {
    hints: false,
  },

  entry: [development && 'webpack/hot/poll?100', './src/index.ts'].filter(Boolean),

  output: {
    filename: outputFilename,
    path: path.join(buildPath, 'node'),
    pathinfo: !development,
  },

  module: {
    rules: [
      ...webpackLoadJavascript({
        babelConfig,
        sourceMaps,
      }),
    ],
  },

  resolve: {
    symlinks: false,
    extensions: ['.wasm', '.ts', '.tsx', '.mjs', '.es.js', '.js', '.jsx', '.json', '.yml'],
    alias: alias(development),
  },

  externals: [
    nodeExternals({
      // @ts-ignore
      allowlist: [development && 'webpack/hot/poll?100'].filter(Boolean),
    }),
  ],

  node: false,

  plugins: [
    new ExtraWatchWebpackPlugin({
      files: [path.join(moduleRoot, 'src/types/**/*.d.ts')],
      dirs: [],
    }),

    webpackTsChecker({
      memoryLimit: 1024,
      eslint: true,
      typeChecking: true,
      exclude: ['src/**/__tests__/**/*.{js,jsx,ts,tsx}', 'src/**/*.(spec|test).{js,jsx,ts,tsx}'],
    }),

    new webpack.EnvironmentPlugin({
      BABEL_ENV: process.env.BABEL_ENV,
      DEBUGGABLE_PROD: process.env.DEBUGGABLE_PROD,
      NODE_ENV: process.env.NODE_ENV,
    }),

    development && new webpack.HotModuleReplacementPlugin(),

    development && new StartServerPlugin({ name: outputFilename }),

    new webpack.optimize.LimitChunkCountPlugin({ maxChunks: 1 }),
  ].filter(Boolean),

  optimization: {
    nodeEnv: false,
    concatenateModules: false,
    noEmitOnErrors: true,
    occurrenceOrder: false,
    removeAvailableModules: true,
    removeEmptyChunks: true,
    minimize: true,
    runtimeChunk: false,
    splitChunks: false,
    minimizer: [production && webpackTerser({ sourceMaps, node: true, profile: false })].filter(Boolean),
  },
}
