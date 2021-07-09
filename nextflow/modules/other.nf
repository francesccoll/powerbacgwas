#!/usr/bin/env nextflow

process copy_files {

    input:
        val(tag)
        each file(file_to_copy)
        path(target_dir)

    """
    cp "${file_to_copy}" "$target_dir"/"${file_to_copy}"
    """
}
