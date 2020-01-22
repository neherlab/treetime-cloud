import path from 'path'

import React, { useCallback } from 'react'
import { connect } from 'react-redux'

import { FileWithPath, useDropzone } from 'react-dropzone'

import { State } from '../../state/reducer'
import { addFiles, AddFilesPayload } from '../../state/upload/upload.actions'
import { selectFiles } from '../../state/upload/upload.selectors'
import { FileType } from '../../state/upload/upload.types'

/* Adds relevant files to a Map to be dispatched */
function reduceDroppedFiles(files: Map<FileType, string>, file: FileWithPath) {
  const type = fileExtToType(path.extname(file.name))
  if (type) {
    files.set(type, file.name)
  }
  return files
}

/* Converts file extension to file FileType enum */
function fileExtToType(ext: string) {
  const extMap = new Map<string, FileType>(
    Object.entries({
      '.nwk': FileType.NWK,
      '.fasta': FileType.FASTA,
      '.csv': FileType.DATES,
      '.json': FileType.CONFIG,
    }),
  )
  return extMap.get(ext)
}

export interface FileUploadZoneProps {
  files: Map<FileType, string>
  addFiles(payload: AddFilesPayload): void
}

function FileUploadZone({ addFiles, files }: FileUploadZoneProps) {
  const onDrop = useCallback(
    (droppedFiles: FileWithPath[]) => {
      const files = droppedFiles.reduce(reduceDroppedFiles, new Map())
      addFiles({ files })
    },
    [addFiles],
  )

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ onDrop })

  return (
    <div {...getRootProps()}>
      <input type="file" {...getInputProps()} />
      {isDragActive ? (
        <p>{'Drop the files here ...'}</p>
      ) : (
        <p>{`Drag 'n' drop some files here, or click to select files`}</p>
      )}
      <ul>
        {[...files.values()].map(filename => (
          <li key={filename}>{filename}</li>
        ))}
      </ul>
    </div>
  )
}

const mapStateToProps = (state: State) => ({
  files: selectFiles(state),
})

const mapDispatchToProps = { addFiles }

export default connect(mapStateToProps, mapDispatchToProps)(FileUploadZone)
