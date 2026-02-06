use anndata::ArrayData;
use anndata::data::DynCsrMatrix;
use anndata::{AnnData, AnnDataOp, ArrayElemOp};
use anndata_hdf5::H5;
use nalgebra_sparse::CsrMatrix;
use num_traits::Zero;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelIterator;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;

pub fn hash_cells(
    adata: AnnData<H5>,
    intersection: Vec<String>,
    hash_seed: i64,
    output_file: &str,
) -> anyhow::Result<()> {
    let mut mask: Vec<usize> = Vec::with_capacity(intersection.len());
    let intersection_set: HashSet<String> = HashSet::from_iter(intersection.into_iter());

    for (idx, name) in adata.var_names().into_iter().enumerate() {
        if intersection_set.contains(&name) {
            mask.push(idx);
        }
    }

    let x = adata.x();
    let n_obs = adata.n_obs();
    let obs_names: Vec<String> = adata.obs_names().into_iter().collect::<Vec<String>>();
    let mut hashes: Vec<(String, [u8; 32])> = Vec::with_capacity(n_obs);

    for (m, start, end) in x.iter::<ArrayData>(512) {
        match m {
            ArrayData::CsrMatrix(dyn_csr) => match dyn_csr {
                DynCsrMatrix::I8(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes)
                }
                DynCsrMatrix::I16(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes)
                }
                DynCsrMatrix::I32(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::I64(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::U8(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::U16(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::U32(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::U64(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::F32(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                DynCsrMatrix::F64(csr) => {
                    hash_csr_range(csr, start, end, &mask, hash_seed, &obs_names, &mut hashes);
                }
                _ => todo!("hashing not implemented for non-numeric type"),
            },
            _ => todo!("This type isnt implemented yet."),
        };
    }

    let mut file = File::options().create(true).write(true).open(output_file)?;
    for (name, hash) in hashes.into_iter() {
        writeln!(&mut file, "{},{}", name, hex::encode(hash))?;
    }
    Ok(())
}

fn hash_csr_range<T: Zero + Sized + Sync + Clone>(
    csr: CsrMatrix<T>,
    start: usize,
    _end: usize,
    mask: &[usize],
    _hash_seed: i64,
    obs_names: &[String],
    ref_buf: &mut Vec<(String, [u8; 32])>,
) {
    let nrows = csr.nrows();

    let name_hash_vec: Vec<(String, [u8; 32])> = (0..nrows)
        .into_par_iter()
        .map(|idx| {
            let g_id = idx + start;
            let hash = csr.hash_row_masked(idx, mask).unwrap();
            let name = obs_names[g_id].clone();
            (name, hash)
        })
        .collect();
    for v in name_hash_vec.into_iter() {
        ref_buf.push(v);
    }
}

fn get_mask_from_intersection(
    gene_names: &[String],
    selected: Vec<String>,
) -> anyhow::Result<Vec<bool>> {
    let sel_set: HashSet<String> = HashSet::from_iter(selected);
    let mut b_vec = vec![false; gene_names.len()];
    for (i, n) in gene_names.iter().enumerate() {
        if sel_set.contains(n) {
            b_vec[i] = true
        }
    }
    Ok(b_vec)
}

trait HashRow {
    fn hash_row(&self, i: usize) -> anyhow::Result<[u8; 32]>;
    fn hash_row_masked(&self, i: usize, mask: &[usize]) -> anyhow::Result<[u8; 32]>;
}

impl<T: Zero + Sized + Clone> HashRow for CsrMatrix<T> {
    fn hash_row(&self, i: usize) -> anyhow::Result<[u8; 32]> {
        let row = self.get_row(i).unwrap();
        let mut d: Vec<T> = Vec::with_capacity(self.ncols());
        for idx in 0..self.ncols() {
            let val = if let Some(v) = row.get_entry(idx) {
                v.into_value()
            } else {
                T::zero()
            };
            d.push(val);
        }
        let byte_slice: &[u8] = unsafe {
            std::slice::from_raw_parts(d.as_ptr() as *const u8, d.len() * size_of::<T>())
        };

        let hash: [u8; 32] = Sha256::digest(byte_slice).into();
        Ok(hash)
    }

    fn hash_row_masked(&self, i: usize, mask: &[usize]) -> anyhow::Result<[u8; 32]> {
        let row = self.get_row(i).unwrap();
        let mut d: Vec<T> = Vec::with_capacity(mask.len());
        for &idx in mask {
            let v = if let Some(x) = row.get_entry(idx) {
                x.into_value()
            } else {
                T::zero()
            };
            d.push(v);
        }
        let byte_slice: &[u8] = unsafe {
            std::slice::from_raw_parts(d.as_ptr() as *const u8, d.len() * size_of::<T>())
        };

        let hash: [u8; 32] = Sha256::digest(byte_slice).into();
        Ok(hash)
    }
}
