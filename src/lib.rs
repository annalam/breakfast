use std::path::Path;
use std::fs::File;
use std::process::Command;

pub fn mkdir(path: &String) {
  if Path::new(path).exists(){
     println!("{:?} already exists", path);
  } else {
    Command::new("mkdir")
            .arg(path)
            .spawn()
            .expect("failed to mkdir");
  }
}

pub fn zopen(path: &String, mode: String) {
  

}
