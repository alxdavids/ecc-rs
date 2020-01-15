#![allow(
    missing_copy_implementations,
    missing_debug_implementations,
    non_camel_case_types,
    non_snake_case,
    unsafe_code
)]
// `#[derive(...)]` uses `trivial_numeric_casts` and `unused_qualifications`
// internally.
#![deny(
    unstable_features, // Used by `internal_benches`
    unused_qualifications,
    variant_size_differences,
)]
// #![forbid(
//     anonymous_parameters,
//     trivial_casts,
//     trivial_numeric_casts,
//     unused_extern_crates,
//     unused_import_braces,
//     unused_results,
//     warnings
// )]
#![no_std]
#![cfg_attr(feature = "internal_benches", allow(unstable_features), feature(test))]

#[cfg(feature = "alloc")]
extern crate alloc;

#[macro_use]
pub mod debug;

#[macro_use]
pub mod test;

#[macro_use]
pub mod arithmetic;

#[macro_use]
mod polyfill;

pub mod bits;
pub(crate) mod c;
pub mod io;
pub mod cpu;

pub mod digest;
pub mod ec;
pub mod endian;
pub mod error;
pub mod limb;
pub mod rand;

mod sealed {
    /// Traits that are designed to only be implemented internally in *ring*.
    //
    // Usage:
    // ```
    // use crate::sealed;
    //
    // pub trait MyType: sealed::Sealed {
    //     // [...]
    // }
    //
    // impl sealed::Sealed for MyType {}
    // ```
    pub trait Sealed {}
}
