declare module '*.json' {
  const content: Record<string, unknown> | Array<unknown>
  export default content
}
