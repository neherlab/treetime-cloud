import React from 'react'
import { connect } from 'react-redux'

import FileUploadZone from '../components/FileUpload/FileUploadZone'

import { State } from '../state/reducer'
import { selectTaskId } from '../state/task/task.selectors'
import { TaskId } from '../state/task/task.types'

import './Home.scss'

export interface HomeProps {
  readonly taskId?: TaskId
}

function Home({ taskId }: HomeProps) {
  return (
    <>
      <h1 className="h1-home">{'Upload files'}</h1>

      <p>{`TaskID: ${taskId}`}</p>

      <div>
        <FileUploadZone />
      </div>
    </>
  )
}

const mapStateToProps = (state: State) => ({
  taskId: selectTaskId(state),
})

export default connect(mapStateToProps)(Home)
