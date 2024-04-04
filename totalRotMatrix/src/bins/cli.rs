#[cfg(not(feature = "cli"))]
compile_error!("Enable --features cli");

use clap::Parser;

#[allow(non_snake_case)]
#[derive(Debug, Parser, serde::Deserialize)]
pub struct Args {
    /// Provide in the form of
    /// `'[1, [2, N], [V00, V01, V10, V11, ..]]'` where `N` is the number of sub-arrays
    /// and `VXY` is the value at the index `[X, Y]`.
    #[clap(short, long, value_parser = SerdeParser::<ndarray::Array2<f64>>::default())]
    PS: ndarray::Array2<f64>,
    #[clap(short, long)]
    T: f64,
    #[clap(short, long)]
    Umax: f64,
    #[clap(short, long)]
    offset: f64,
    #[clap(short, long)]
    inpoFact: u64,
    /// Provide in the form of
    /// `'[1, [N], [V0, V1, ..]]'` where `N` is the length of the array and `VX` is the
    /// value at index `[X]`
}

fn main() -> serde_json::Result<()> {
    let args = Args::parse();
    let Rtot = totalRotMatrix::create_Rtot(
        args.PS.view(),
        args.T,
        args.Umax,
        args.offset,
        args.inpoFact,
    );

    serde_json::to_writer(std::io::stdout(), &Rtot)
}

#[derive(Clone, Default)]
struct SerdeParser<T>(std::marker::PhantomData<T>);

impl<T> clap::builder::TypedValueParser for SerdeParser<T>
where
    T: Clone + for<'de> serde::Deserialize<'de> + Send + Sync + 'static,
{
    type Value = T;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<Self::Value, clap::Error> {
        let s = clap::builder::NonEmptyStringValueParser::new().parse_ref(cmd, arg, value)?;
        Ok(serde_json::from_str(&s).unwrap())
    }
}
