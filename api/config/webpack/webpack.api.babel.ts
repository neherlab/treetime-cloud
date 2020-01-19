import '../dotenv'

import path from 'path'

import ExtraWatchWebpackPlugin from 'extra-watch-webpack-plugin'
import NodeHotLoaderWebpackPlugin from 'node-hot-loader/NodeHotLoaderWebpackPlugin'
import webpack from 'webpack'
import nodeExternals from 'webpack-node-externals'

import webpackFriendlyConsole from './lib/webpackFriendlyConsole'
import webpackLoadJavascript from './lib/webpackLoadJavascript'
import webpackTerser from './lib/webpackTerser'
import webpackTsChecker from './lib/webpackTsChecker'

import { findModuleRoot } from '../../lib/findModuleRoot'
import { getenv } from '../../src/lib/getenv'

import babelConfig from '../../babel.config'

const MODE: 'development' | 'production' = getenv('NODE_ENV') === 'development' ? 'development' : 'production' // prettier-ignore

const production = MODE === 'production'
const development = MODE === 'development'
const sourceMaps = true
const fancyConsole = getenv('DEV_FANCY_CONSOLE', '0') === '1'
const fancyClearConsole = getenv('DEV_FANCY_CLEAR_CONSOLE', '0') === '1'

const { moduleRoot, pkg } = findModuleRoot()
const buildPath = path.join(moduleRoot, '.build', MODE)

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
  stats: fancyConsole
    ? false
    : {
        all: false,
        errors: true,
        warnings: true,
        moduleTrace: true,
        colors: true,
      },
  performance: {
    hints: false,
  },

  entry: `./src/index.${MODE}`,

  output: {
    filename: 'api.js',
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
    extensions: [
      '.wasm',
      '.ts',
      '.tsx',
      '.mjs',
      '.es.js',
      '.js',
      '.jsx',
      '.json',
      '.yml',
    ],
    alias: alias(development),
  },

  externals: [nodeExternals()],

  node: {
    __dirname: true,
  },

  plugins: [
    new ExtraWatchWebpackPlugin({
      files: [path.join(moduleRoot, 'src/types/**/*.d.ts')],
      dirs: [],
    }),

    webpackTsChecker({
      memoryLimit: 1024,
      tslint: path.join(moduleRoot, 'tslint.json'),
      tsconfig: path.join(moduleRoot, 'tsconfig.json'),
      reportFiles: [
        'src/**/*.{js,jsx,ts,tsx}',
        '!src/**/__tests__/**/*.{js,jsx,ts,tsx}',
        '!src/**/*.(spec|test).{js,jsx,ts,tsx}',
        '!static/**/*',
      ],
    }),

    new webpack.EnvironmentPlugin({
      BABEL_ENV: process.env.BABEL_ENV,
      DEBUGGABLE_PROD: process.env.DEBUGGABLE_PROD,
      NODE_ENV: process.env.NODE_ENV,
    }),

    ...(fancyConsole
      ? webpackFriendlyConsole({
          clearConsole: fancyClearConsole,
          projectRoot: path.resolve(moduleRoot, '..'),
          packageName: pkg.name || 'api',
          progressBarColor: 'red',
        })
      : []),

    development && new NodeHotLoaderWebpackPlugin({force: true, logLevel: 'minimal'}), // prettier-ignore

    development && new webpack.HotModuleReplacementPlugin(),

    new webpack.optimize.LimitChunkCountPlugin({ maxChunks: 1 }),
  ].filter(Boolean),

  optimization: {
    concatenateModules: false,
    noEmitOnErrors: true,
    occurrenceOrder: false,
    removeAvailableModules: true,
    removeEmptyChunks: true,
    minimize: true,
    runtimeChunk: false,
    splitChunks: false,
    minimizer: [
      production && webpackTerser({ sourceMaps, node: true, profile: false }),
    ].filter(Boolean),
  },
}
