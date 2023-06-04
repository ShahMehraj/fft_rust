/**
 * @file main.c
 * @brief WAV file processing using FFT
 *
 * This program reads audio samples from a WAV file, performs Fast Fourier Transform (FFT) on the samples,
 * and prints the real and imaginary parts of the FFT result.
 *
 * The program processes the audio samples in blocks of size `BLOCK_SIZE`. Processing in larger block sizes
 * may lead to unresponsiveness, so it is recommended to keep `BLOCK_SIZE` as 1024 or lower.
 *
 * For the remaining samples that are not a multiple of `BLOCK_SIZE`, the program finds the largest power of 4
 * that is less than or equal to the remaining samples and performs FFT on that size.
 *
 *
 * Note: This implementation assumes the availability of the `read_wav` function from the "wave.h" library,
 * as well as the `impeghd_rad2_cplx_fft` function from the "fft.h" library. Make sure to include the correct
 * header files and provide the necessary library files when compiling the program.
 *
 * @date 2023-06-04
 **/
mod fft;
mod read_wav;
use fft::impeghd_rad2_cplx_fft;
use read_wav::read_wav;

fn main() {
    let file_path = "1khz_Sine_44_1khz.wav";
    let (num_samples, mut ptr_real) = read_wav(file_path);
    //println!("{:?}",ptr_real);

    let block_size = 1024;
    let num_of_iterations = num_samples / block_size;
    //num_samples = num_of_iterations * block_size;
    let mut ptr_imag: Vec<f32> = vec![0.0; num_samples as usize];
    let mut ptr_scratch: Vec<f32> = vec![0.0; (4 * 1024) as usize];

    for i in 0..num_of_iterations {
        let start_index = (i * block_size) as usize;
        let end_index = start_index + block_size as usize;

        impeghd_rad2_cplx_fft(
            &mut ptr_real[start_index..end_index],
            &mut ptr_imag[start_index..end_index],
            block_size as u32,
            &mut ptr_scratch,
        );
    }

    /*
     The following block finds the largest
     power of 4 less than or equal to the
     remaining_samples and computes their fft
    */

    let remaining_samples = num_samples % block_size;
    if remaining_samples > 0 {
        let start_index = (num_samples - remaining_samples) as usize;
    
        let mut power_of_4 = 1;
        while power_of_4 * 4 <= remaining_samples {
            power_of_4 *= 4;
            
        }
        
        let end_index = start_index + power_of_4 as usize;

        impeghd_rad2_cplx_fft(
            &mut ptr_real[start_index..end_index],
            &mut ptr_imag[start_index..end_index],
            power_of_4 as u32,
            &mut ptr_scratch,
        );
    }
    println!("SUCESSES!");
    for i in 0..num_samples
    {
        let c = if ptr_imag[i as usize] < 0.0 { '-' } else { '+' };
        println!("{} {} {} i", ptr_real[i as usize], c , ptr_imag[i as usize].abs());
    }
}

