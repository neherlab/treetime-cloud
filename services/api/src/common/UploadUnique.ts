import { UseInterceptors } from '@nestjs/common'
import { FileFieldsInterceptor } from '@nestjs/platform-express'

/**
 * Decorator. Uses `FileFieldsInterceptor` to verify that form data contains
 * one file upload field of a given type at most.
 *
 * @param fields Fields in the form data to be verified
 */
export function UploadUnique(fields: string[]) {
  return UseInterceptors(
    FileFieldsInterceptor(
      fields.map((field) => ({ name: field, maxCount: 1 })),
    ),
  )
}
