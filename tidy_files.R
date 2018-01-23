source("R/utilities.R")


cmd_arg <- commandArgs(trailingOnly = T)
sex <- cmd_arg[1]
config_file <- cmd_arg[2]
config <- yaml::read_yaml(config_file)

# copy files to a time-stamped folder inside the existing results folder.
# (for archiving purposes)
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

base_dir <- config$out_dir

out_path <- file.path(base_dir, config$type, tolower(sex))
archive_path <- file.path(out_path, time_stamp)
dir.create(archive_path)

file_list <- get_non_directories(out_path)

purrr::walk(file_list,
            function(file, out_path, archive_path){
              file.copy(file.path(out_path, file), archive_path)
              # file.remove(file.path(out_path, file))
            },
            out_path, archive_path)

# in case you wish to add an additional script to copy outputs elsewhere.
if (file.exists("backup_copy.sh")){
  command <- paste0("backup_copy.sh ", archive_path)
}