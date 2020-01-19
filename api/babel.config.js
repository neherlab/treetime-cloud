require('./config/dotenv')

const loose = true

module.exports = {
  compact: false,
  presets: [
    '@babel/preset-typescript',
    [
      '@babel/preset-env',
      {
        loose,
        corejs: false,
        modules: 'commonjs',
        shippedProposals: true,
        targets: { node: '12' },
        exclude: ['transform-typeof-symbol'],
      },
    ],
  ],
  plugins: [
    'babel-plugin-transform-typescript-metadata', // goes before "proposal-decorators"
    ['@babel/plugin-proposal-decorators', { legacy: true }], // goes before "class-properties"
    'babel-plugin-parameter-decorator',
    ['@babel/plugin-proposal-class-properties', { loose }],
    ['@babel/plugin-proposal-numeric-separator', { loose }],
  ].filter(Boolean),
}
