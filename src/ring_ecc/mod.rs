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
    // use crate::ring_ecc::sealed;
    //
    // pub trait MyType: sealed::Sealed {
    //     // [...]
    // }
    //
    // impl sealed::Sealed for MyType {}
    // ```
    pub trait Sealed {}
}
