use std::path::Path;
use anyhow::Result;

/// Types implementing this trait can be saved to file.
pub trait Save {
    /// Serialise the type to a given file
    /// # Errors
    /// if the instance can not be serialised or if the file can't be written to.
    fn save_data(&self, path: &Path) -> Result<()>;

    /// Report the saving of a file (if it is a filepath) and save the data.
    /// # Errors
    /// if the instance can not be serialised or if the file can't be written to.
    fn save(&self, path: &Path) -> Result<()> {
        println!("[SAVE] {}", path.display());

        self.save_data(path)
    }
}
