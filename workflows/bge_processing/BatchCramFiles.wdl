version 1.0

# Coerces one manifest slice to Array[File] so Cromwell localizes only this batch's paths.
workflow BatchCramFiles {
  input {
    Array[Array[String]] manifest_rows
    Int header_offset
    Int batch_idx
    Int batch_size
    Int total_crams
  }

  Int row_start = batch_idx * batch_size
  Int n_this = if (batch_idx + 1) * batch_size < total_crams then batch_size else total_crams - row_start

  scatter (i in range(n_this)) {
    File cram = manifest_rows[header_offset + row_start + i][0]
    File crai = manifest_rows[header_offset + row_start + i][1]
  }

  output {
    Array[File] crams = cram
    Array[File] crais = crai
  }
}
