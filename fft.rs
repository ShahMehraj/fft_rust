pub fn ia_negate_flt(a: f32) -> f32 {
    -a
}

pub fn ia_sub_flt(a: f32, b: f32) -> f32 {
    a - b
}

pub fn ia_add_flt(a: f32, b: f32) -> f32 {
    a + b
}

pub fn ia_mul_flt(a: f32, b: f32) -> f32 {
    a * b
}

pub fn ia_mac_flt(x: f32, a: f32, b: f32) -> f32 {
    x + a * b
}

pub fn ia_msu_flt(x: f32, a: f32, b: f32) -> f32 {
    x - a * b
}

pub fn norm32(a: u32) -> u32 {
    let mut norm_val: u32;

    if a == 0 {
        norm_val = 31;
    } else if a == 0xffffffff {
        norm_val = 31;
    } else {
        let mut a = a as i32;
        if a < 0 {
            a = !a;
        }
        norm_val = 0;
        while a < 0x40000000 {
            a <<= 1;
            norm_val += 1;
        }
    }

    norm_val
}
pub fn dig_rev(i: u32, m: u32) -> u32 {
    let mut result = i;
    result = ((result & 0x33333333) << 2) | ((result & !0x33333333) >> 2);
    result = ((result & 0x0F0F0F0F) << 4) | ((result & !0x0F0F0F0F) >> 4);
    result = ((result & 0x00FF00FF) << 8) | ((result & !0x00FF00FF) >> 8);
    result >> m
}
pub fn impeghd_rad2_cplx_fft(
    ptr_real: &mut [f32],
    ptr_imag: &mut [f32],
    n_points: u32,
    ptr_scratch: &mut [f32],
) {
    let mut j: u32;
    let mut n_stages: u32;
    let mut x0r: f32;
    let mut x0i: f32;
    let mut x1r: f32;
    let mut x1i: f32;
    let mut x2r: f32;
    let mut x2i: f32;
    let mut x3r: f32;
    let mut x3i: f32;
    let mut del: u32;
    let mut nodespacing: u32;
    let mut in_loop_cnt: u32;
    let not_power_4: bool;
    let dig_rev_shift: u32;
    let m_points = n_points;
    let ptr_x: *mut f32 = ptr_scratch.as_mut_ptr();
    let y: *mut f32 = unsafe { ptr_scratch.as_mut_ptr().add(2048) };
    let mut ptr_y = y;    

    dig_rev_shift = norm32(m_points as u32) + 1 - 16;
    n_stages = 30 - norm32(m_points as u32);
    not_power_4 = n_stages & 1 != 0;
    n_stages = n_stages >> 1;

    let ptr_w: *const f32 = IA_FFT_TWIDDLE_TABLE_FLOAT.as_ptr();

    for i in 0..n_points {
        unsafe {
            *ptr_x.offset(2 * i as isize) = ptr_real[i as usize];
            *ptr_x.offset(2 * i as isize + 1) = ptr_imag[i as usize];
        }
    }
       
    

    for i in (0..n_points).step_by(4)
    {
      let mut tmk : f32;
      let mut inp = ptr_x;
  
      let mut h2 = dig_rev(i as u32, dig_rev_shift);
      if not_power_4 {
          h2 += 1;
          h2 &= !1;
      }
      inp = unsafe{ inp.add(h2 as usize) };
  
      x0r = unsafe { *inp }; //println!("{}", x0r);
      x0i = unsafe { *inp.offset(1) };
      inp = unsafe{ inp.add((n_points >> 1) as usize) };

      x1r = unsafe{ *inp }; //println!("{}", x1r);
      x1i = unsafe { *inp.offset(1) };
      inp = unsafe{ inp.add((n_points >> 1) as usize) };

      x2r = unsafe{ *inp }; //println!("{}", x2r);
      x2i = unsafe { *inp.offset(1) };
      inp = unsafe{ inp.add((n_points >> 1) as usize) };

      x3r = unsafe{ *inp }; //println!("{}", x3r);
      x3i = unsafe { *inp.offset(1) };
  
      x0r = ia_add_flt(x0r, x2r);
      x0i = ia_add_flt(x0i, x2i);
  
      tmk = ia_sub_flt(x0r, x2r);
      x2r = ia_sub_flt(tmk, x2r);
      tmk = ia_sub_flt(x0i, x2i);
      x2i = ia_sub_flt(tmk, x2i);
  
      x1r = ia_add_flt(x1r, x3r);
      x1i = ia_add_flt(x1i, x3i);
  
      tmk = ia_sub_flt(x1r, x3r);
      x3r = ia_sub_flt(tmk, x3r);
      tmk = ia_sub_flt(x1i, x3i);
      x3i = ia_sub_flt(tmk, x3i);
  
      x0r = ia_add_flt(x0r, x1r);
      x0i = ia_add_flt(x0i, x1i);
  
      tmk = ia_sub_flt(x0r, x1r);
      x1r = ia_sub_flt(tmk, x1r);
      tmk = ia_sub_flt(x0i, x1i);
      x1i = ia_sub_flt(tmk, x1i);
  
      x2r = ia_add_flt(x2r, x3i);
      x2i = ia_sub_flt(x2i, x3r);
  
      tmk = ia_sub_flt(x2r, x3i);
      x3i = ia_sub_flt(tmk, x3i);
      tmk = ia_add_flt(x2i, x3r);
      x3r = ia_add_flt(tmk, x3r);
  
      unsafe
      {
          *ptr_y = x0r;
          ptr_y = ptr_y.add(1);
          *ptr_y = x0i;
          ptr_y = ptr_y.add(1);
          *ptr_y = x2r;
          ptr_y = ptr_y.add(1);
          *ptr_y = x2i;
          ptr_y = ptr_y.add(1);
          *ptr_y = x1r;
          ptr_y = ptr_y.add(1);
          *ptr_y = x1i;
          ptr_y = ptr_y.add(1);
          *ptr_y = x3i;
          ptr_y = ptr_y.add(1);
          *ptr_y = x3r;
          ptr_y = ptr_y.add(1);
      }
       
    }
    
    unsafe {
        //ptr_y = ptr_y.sub(2 * n_points as isize);
        ptr_y = ptr_y.offset(-(2 * n_points as isize));
    }

    del = 4;
    nodespacing = 64;
    in_loop_cnt = (n_points >> 4) as u32;
    
    for _i in (1..n_stages).rev() {

        let twiddles: *const f32 = ptr_w;
        let mut data: *mut f32 = ptr_y;
        let (mut w1, mut w2, mut w3, mut w4, mut w5, mut w6): (f32, f32, f32, f32, f32, f32);
        let mut sec_loop_cnt: u32;

        for _k in (1..=in_loop_cnt).rev() {
            x0r = unsafe { *data };
            x0i = unsafe { *(data.add(1)) };
            data = unsafe { data.add((del << 1) as usize) };

            x1r = unsafe { *data };
            x1i = unsafe { *(data.add(1)) };
            data = unsafe { data.add((del << 1) as usize) };

            x2r = unsafe { *data };
            x2i = unsafe { *(data.add(1)) };
            data = unsafe { data.add((del << 1) as usize) };

            x3r = unsafe { *data };
            x3i = unsafe { *(data.add(1)) };
            data = unsafe { data.sub((3 * (del << 1)) as usize) };

            x0r = ia_add_flt(x0r, x2r);
            x0i = ia_add_flt(x0i, x2i);
            x2r = ia_msu_flt(x0r, x2r, 2.0);
            x2i = ia_msu_flt(x0i, x2i, 2.0);
            x1r = ia_add_flt(x1r, x3r);
            x1i = ia_add_flt(x1i, x3i);
            x3r = ia_msu_flt(x1r, x3r, 2.0);
            x3i = ia_msu_flt(x1i, x3i, 2.0);
      
            x0r = ia_add_flt(x0r, x1r);
            x0i = ia_add_flt(x0i, x1i);
            x1r = ia_msu_flt(x0r, x1r, 2.0);
            x1i = ia_msu_flt(x0i, x1i, 2.0);
            x2r = ia_add_flt(x2r, x3i);
            x2i = ia_sub_flt(x2i, x3r);
            x3i = ia_msu_flt(x2r, x3i, 2.0);
            x3r = ia_mac_flt(x2i, x3r, 2.0);

            unsafe {
                *data = x0r;
                *(data.add(1)) = x0i;
                data = data.add((del << 1) as usize);

                *data = x2r;
                *(data.add(1)) = x2i;
                data = data.add((del << 1) as usize);

                *data = x1r;
                *(data.add(1)) = x1i;
                data = data.add((del << 1) as usize);

                *data = x3i;
                *(data.add(1)) = x3r;
                data = data.add((del << 1) as usize);
            }
        }
        unsafe{data = ptr_y.add(2);}
        sec_loop_cnt = nodespacing * del;
        sec_loop_cnt = (sec_loop_cnt / 4) + (sec_loop_cnt / 8) - (sec_loop_cnt / 16)
            + (sec_loop_cnt / 32)
            - (sec_loop_cnt / 64)
            + (sec_loop_cnt / 128)
            - (sec_loop_cnt / 256);
        j = nodespacing;
        //for j in (nodespacing..=sec_loop_cnt).step_by(nodespacing as usize) 
        while j <= sec_loop_cnt
        {
            
            w1 = unsafe{*twiddles.offset(j as isize)};
            w4 = unsafe{*twiddles.offset(j as isize + 257)};
            w2 = unsafe{*twiddles.offset((j << 1) as isize)};
            w5 = unsafe{*twiddles.offset((j << 1) as isize + 257)};
            w3 = unsafe{*twiddles.offset((j + (j << 1)) as isize)};
            w6 = unsafe{*twiddles.offset((j + (j << 1)) as isize + 257)};
            

            // twiddles += nodespacing;

            for _k in (1..=in_loop_cnt).rev() {
                
				let mut tmp:  f32;
				unsafe {
                    
                    // x0 is loaded later to avoid register crunch

                    data = data.offset((del << 1) as isize);

                    x1r = *data;
                    x1i = *(data.offset(1));
                    data = data.offset((del << 1) as isize);

                    x2r = *data;
                    x2i = *(data.offset(1));
                    data = data.offset((del << 1) as isize);

                    x3r = *data;
                    x3i = *(data.offset(1));
                    let offset = (3 * (del << 1)) as usize;
                    data = data.offset(offset.wrapping_neg() as isize);

                }

                tmp = ia_sub_flt(ia_mul_flt(x1r, w1), ia_mul_flt(x1i, w4));
                x1i = ia_mac_flt(ia_mul_flt(x1r, w4), x1i, w1);
                x1r = tmp;
        
                tmp = ia_sub_flt(ia_mul_flt(x2r, w2), ia_mul_flt(x2i, w5));
                x2i = ia_mac_flt(ia_mul_flt(x2r, w5), x2i, w2);
                x2r = tmp;
        
                tmp = ia_sub_flt(ia_mul_flt(x3r, w3), ia_mul_flt(x3i, w6));
                x3i = ia_mac_flt(ia_mul_flt(x3r, w6), x3i, w3);
                x3r = tmp;

                unsafe
                {
                    x0r = *data;
                    x0i = *data.offset(1); 
                }

                x0r = ia_add_flt(x0r, x2r);
                x0i = ia_add_flt(x0i, x2i);
                x2r = ia_msu_flt(x0r, x2r, 2.0);
                x2i = ia_msu_flt(x0i, x2i, 2.0);
                x1r = ia_add_flt(x1r, x3r);
                x1i = ia_add_flt(x1i, x3i);
                x3r = ia_msu_flt(x1r, x3r, 2.0);
                x3i = ia_msu_flt(x1i, x3i, 2.0);
        
                x0r = ia_add_flt(x0r, x1r);
                x0i = ia_add_flt(x0i, x1i);
                x1r = ia_msu_flt(x0r, x1r, 2.0);
                x1i = ia_msu_flt(x0i, x1i, 2.0);
                x2r = ia_add_flt(x2r, x3i);
                x2i = ia_sub_flt(x2i, x3r);
                x3i = ia_msu_flt(x2r, x3i, 2.0);
                x3r = ia_mac_flt(x2i, x3r, 2.0);

                unsafe {
                    *data = x0r;
                    *(data.offset(1)) = x0i;
                    data = data.offset((del << 1) as isize);

                    *data = x2r;
                    *(data.offset(1)) = x2i;
                    data = data.offset((del << 1) as isize);

                    *data = x1r;
                    *(data.offset(1)) = x1i;
                    data = data.offset((del << 1) as isize);

                    *data = x3i;
                    *(data.offset(1)) = x3r;
                    data = data.offset((del << 1) as isize);
                }
            }
            unsafe {
                data = data.offset(-((2 * n_points) as isize)); 
                data = data.offset(2);
            }
            j += nodespacing;
            
            
        }
        //println!("{}", j);
        //let final_j = j;
        //for j in (final_j..=((nodespacing * del) >> 1)).step_by(nodespacing as usize) 
        while j <= (nodespacing * del) >> 1
        {
            w1 = unsafe{*twiddles.offset(j as isize)};
            w4 = unsafe{*twiddles.offset(j as isize + 257)};
            w2 = unsafe{*twiddles.offset((j << 1) as isize)};
            w5 = unsafe{*twiddles.offset((j << 1) as isize + 257)};
            w3 = unsafe{*twiddles.offset((j + (j << 1)) as isize -256)};
            w6 = unsafe{*twiddles.offset((j + (j << 1)) as isize + 1)};

            // twiddles += nodespacing;

            for _k in (1..=in_loop_cnt).rev() {
                let mut tmp:  f32;
                unsafe {
                    data = data.add((del << 1) as usize);

                    x1r = *data;
                    x1i = *data.add(1);
                    data = data.add((del << 1) as usize);

                    x2r = *data;
                    x2i = *data.add(1);
                    data = data.add((del << 1) as usize);

                    x3r = *data;
                    x3i = *data.add(1);
                    data = data.sub((3 * (del << 1)) as usize);
                }

                tmp = ia_sub_flt(ia_mul_flt(x1r, w1), ia_mul_flt(x1i, w4));
                x1i = ia_mac_flt(ia_mul_flt(x1r, w4), x1i, w1);
                x1r = tmp;
        
                tmp = ia_sub_flt(ia_mul_flt(x2r, w2), ia_mul_flt(x2i, w5));
                x2i = ia_mac_flt(ia_mul_flt(x2r, w5), x2i, w2);
                x2r = tmp;
        
                tmp = ia_add_flt(ia_mul_flt(x3r, w6), ia_mul_flt(x3i, w3));
                x3i = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, w3)), ia_mul_flt(x3i, w6));
                x3r = tmp;

                unsafe {
                    x0r = *data;
                    x0i = *data.add(1);
                }
                x0r = ia_add_flt(x0r, x2r);
                x0i = ia_add_flt(x0i, x2i);
                x2r = ia_msu_flt(x0r, x2r, 2.0);
                x2i = ia_msu_flt(x0i, x2i, 2.0);
                x1r = ia_add_flt(x1r, x3r);
                x1i = ia_add_flt(x1i, x3i);
                x3r = ia_msu_flt(x1r, x3r, 2.0);
                x3i = ia_msu_flt(x1i, x3i, 2.0);
        
                x0r = ia_add_flt(x0r, x1r);
                x0i = ia_add_flt(x0i, x1i);
                x1r = ia_msu_flt(x0r, x1r, 2.0);
                x1i = ia_msu_flt(x0i, x1i, 2.0);
                x2r = ia_add_flt(x2r, x3i);
                x2i = ia_sub_flt(x2i, x3r);
                x3i = ia_msu_flt(x2r, x3i, 2.0);
                x3r = ia_mac_flt(x2i, x3r, 2.0);

                unsafe {
                    *data = x0r;
                    *(data.add(1)) = x0i;
                    data = data.add((del << 1) as usize);

                    *data = x2r;
                    *(data.add(1)) = x2i;
                    data = data.add((del << 1) as usize);

                    *data = x1r;
                    *(data.add(1)) = x1i;
                    data = data.add((del << 1) as usize);

                    *data = x3i;
                    *(data.add(1)) = x3r;
                    data = data.add((del << 1) as usize);
                }
                 
            }
            unsafe {
                data = data.sub((2 * n_points) as usize);
                data = data.add(2);
            }
            j += nodespacing;
        }
        //println!("{}", j);
        
        //for j in (((nodespacing * del) >> 1)..=(sec_loop_cnt * 2)).step_by(nodespacing as usize)
        while j <= sec_loop_cnt * 2 
        {
            w1 = unsafe{*twiddles.offset(j as isize)};
            w4 = unsafe{*twiddles.offset(j as isize + 257)};
            w2 = unsafe{*twiddles.offset((j << 1) as isize -256)};
            w5 = unsafe{*twiddles.offset((j << 1) as isize + 1)};
            w3 = unsafe{*twiddles.offset((j + (j << 1)) as isize - 256)};
            w6 = unsafe{*twiddles.offset((j + (j << 1)) as isize + 1)};

            // twiddles += nodespacing;
            for _k in (1..=in_loop_cnt).rev() {
                let mut tmp: f32;
                /*x0 is loaded later to avoid register crunch*/

                data = data.wrapping_add((del << 1) as usize);
                
                x1r = unsafe{*data};
                x1i = unsafe { *(data.add(1)) };
                data = data.wrapping_add((del << 1) as usize);

                x2r = unsafe{*data};
                x2i = unsafe { *(data.add(1)) };
                data = data.wrapping_add((del << 1) as usize);

                x3r = unsafe{*data};
                x3i = unsafe { *(data.add(1)) };
                data = data.wrapping_sub((3 * (del << 1)) as usize);

                tmp = ia_sub_flt(ia_mul_flt(x1r, w1), ia_mul_flt(x1i, w4));
                x1i = ia_mac_flt(ia_mul_flt(x1r, w4), x1i, w1);
                x1r = tmp;
        
                tmp = ia_add_flt(ia_mul_flt(x2r, w5), ia_mul_flt(x2i, w2));
                x2i = ia_add_flt(ia_negate_flt(ia_mul_flt(x2r, w2)), ia_mul_flt(x2i, w5));
                x2r = tmp;
        
                tmp = ia_add_flt(ia_mul_flt(x3r, w6), ia_mul_flt(x3i, w3));
                x3i = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, w3)), ia_mul_flt(x3i, w6));
                x3r = tmp;

                x0r = unsafe { *data };
                x0i = unsafe { *(data.offset(1)) };

                x0r = ia_add_flt(x0r, x2r);
                x0i = ia_add_flt(x0i, x2i);
                x2r = ia_msu_flt(x0r, x2r, 2.0);
                x2i = ia_msu_flt(x0i, x2i, 2.0);
                x1r = ia_add_flt(x1r, x3r);
                x1i = ia_add_flt(x1i, x3i);
                x3r = ia_msu_flt(x1r, x3r, 2.0);
                x3i = ia_msu_flt(x1i, x3i, 2.0);
        
                x0r = ia_add_flt(x0r, x1r);
                x0i = ia_add_flt(x0i, x1i);
                x1r = ia_msu_flt(x0r, x1r, 2.0);
                x1i = ia_msu_flt(x0i, x1i, 2.0);
                x2r = ia_add_flt(x2r, x3i);
                x2i = ia_sub_flt(x2i, x3r);
                x3i = ia_msu_flt(x2r, x3i, 2.0);
                x3r = ia_mac_flt(x2i, x3r, 2.0);

                unsafe {
                    *data = x0r;
                    *(data.offset(1)) = x0i;
                    data = data.offset((del << 1) as isize);

                    *data = x2r;
                    *(data.offset(1)) = x2i;
                    data = data.offset((del << 1) as isize);

                    *data = x1r;
                    *(data.offset(1)) = x1i;
                    data = data.offset((del << 1) as isize);

                    *data = x3i;
                    *(data.offset(1)) = x3r;
                    data = data.offset((del << 1) as isize);
                }
            }
            unsafe {
                data = data.sub((2 * n_points) as usize);
                data = data.add(2);
            }
            j += nodespacing
        }
        //println!("{}", j);
        
        //for j in ((sec_loop_cnt * 2)..=(nodespacing * del)).step_by(nodespacing as usize) 
        while j < nodespacing * del
        {
            w1 = unsafe{*twiddles.offset(j as isize)};
            w4 = unsafe{*twiddles.offset(j as isize + 257)};
            w2 = unsafe{*twiddles.offset((j << 1) as isize -256)};
            w5 = unsafe{*twiddles.offset((j << 1) as isize + 1)};
            w3 = unsafe{*twiddles.offset((j + (j << 1)) as isize - 512)};
            w6 = unsafe{*twiddles.offset((j + (j << 1)) as isize - 512 + 257)};

            // twiddles += nodespacing;
            for _k in (1..=in_loop_cnt).rev() {
                let mut tmp: f32;
                /*x0 is loaded later to avoid register crunch*/

                data = unsafe { data.offset((del * 2) as isize) }; 

                x1r = unsafe { *data };
                x1i = unsafe { *(data.add(1)) };
                data = unsafe { data.offset((del * 2) as isize) };

                x2r = unsafe { *data };
                x2i = unsafe { *(data.add(1)) };
                data = unsafe { data.offset((del * 2) as isize) };

                x3r = unsafe { *data };
                x3i = unsafe { *(data.add(1)) };
                data = unsafe { data.offset(-((del * 6) as isize)) };

                tmp = ia_sub_flt(ia_mul_flt(x1r, w1), ia_mul_flt(x1i, w4));
                x1i = ia_mac_flt(ia_mul_flt(x1r, w4), x1i, w1);
                x1r = tmp;
        
                tmp = ia_add_flt(ia_mul_flt(x2r, w5), ia_mul_flt(x2i, w2));
                x2i = ia_add_flt(ia_negate_flt(ia_mul_flt(x2r, w2)), ia_mul_flt(x2i, w5));
                x2r = tmp;
        
                tmp = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, w3)), ia_mul_flt(x3i, w6));
                x3i = ia_mac_flt(ia_mul_flt(x3r, w6), x3i, w3);
                x3r = tmp;

                x0r = unsafe { *data };
                x0i = unsafe { *(data.add(1)) };

                x0r = ia_add_flt(x0r, x2r);
                x0i = ia_add_flt(x0i, x2i);
                x2r = ia_msu_flt(x0r, x2r, 2.0);
                x2i = ia_msu_flt(x0i, x2i, 2.0);
                x1r = ia_add_flt(x1r, x3r);
                x1i = ia_sub_flt(x1i, x3i);
                x3r = ia_msu_flt(x1r, x3r, 2.0);
                x3i = ia_mac_flt(x1i, x3i, 2.0);
        
                x0r = ia_add_flt(x0r, x1r);
                x0i = ia_add_flt(x0i, x1i);
                x1r = ia_msu_flt(x0r, x1r, 2.0);
                x1i = ia_msu_flt(x0i, x1i, 2.0);
                x2r = ia_add_flt(x2r, x3i);
                x2i = ia_sub_flt(x2i, x3r);
                x3i = ia_msu_flt(x2r, x3i, 2.0);
                x3r = ia_mac_flt(x2i, x3r, 2.0);

                unsafe {
                    *data = x0r;
                    *(data.add(1)) = x0i;
                    data = data.add((del << 1) as usize);

                    *data = x2r;
                    *(data.add(1)) = x2i;
                    data = data.add((del << 1) as usize);

                    *data = x1r;
                    *(data.add(1)) = x1i;
                    data = data.add((del << 1) as usize);

                    *data = x3i;
                    *(data.add(1)) = x3r;
                    data = data.add((del << 1) as usize);
                }
                //println!("{} {} {} {} {} {} {} {}", x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i);
            }
            unsafe {
                data = data.sub((2 * n_points) as usize);
                data = data.add(2);
            }
            j += nodespacing;
        }
        nodespacing >>= 2;
        del <<= 2;
        in_loop_cnt >>= 2;
        
    }


    if not_power_4 {
        let mut twiddles: *mut f32 = ptr_w as *mut f32;
        nodespacing <<= 1;

        for _j in (1..=(del / 2)).rev() {
            let w1 = unsafe { *twiddles };
            let w4 = unsafe { *(twiddles.add(257)) };
            let tmp: f32;

            twiddles = unsafe { twiddles.add(nodespacing as usize) };

            x0r = unsafe { *ptr_y };
            x0i = unsafe { *(ptr_y.add(1)) };
            ptr_y = unsafe { ptr_y.add((del << 1) as usize) };

            x1r = unsafe { *ptr_y };
            x1i = unsafe { *(ptr_y.add(1)) };

            tmp = ia_sub_flt(ia_mul_flt(x1r, w1), ia_mul_flt(x1i, w4));
            x1i = ia_mac_flt(ia_mul_flt(x1r, w4), x1i, w1) as f32;
            x1r = tmp;

            unsafe {
                *ptr_y = ia_sub_flt(x0r, x1r);
                *(ptr_y.add(1)) = ia_sub_flt(x0i, x1i);
                ptr_y = ptr_y.sub((del << 1) as usize); 

                *ptr_y = ia_add_flt(x0r, x1r);
                *(ptr_y.add(1)) = ia_add_flt(x0i, x1i);
                ptr_y = ptr_y.add(2);
            }
        }

        twiddles = ptr_w as *mut f32;
        for _j in (1..=(del / 2)).rev() {
            let w1 = unsafe { *twiddles };
            let w4 = unsafe { *(twiddles.add(257)) };
            let tmp: f32;

            twiddles = unsafe { twiddles.add(nodespacing as usize) };

            x0r = unsafe { *ptr_y };
            x0i = unsafe { *(ptr_y.add(1)) };
            ptr_y = unsafe { ptr_y.add((del << 1) as usize) };

            x1r = unsafe { *ptr_y };
            x1i = unsafe { *(ptr_y.add(1)) };

            tmp = ia_add_flt(ia_mul_flt(x1r, w4), ia_mul_flt(x1i, w1));
            x1i = ia_add_flt(ia_negate_flt(ia_mul_flt(x1r, w1)), ia_mul_flt(x1i, w4));
            x1r = tmp;

            unsafe {
                *ptr_y = ia_sub_flt(x0r, x1r);
                *(ptr_y.add(1)) = ia_sub_flt(x0i, x1i);
                ptr_y = ptr_y.sub((del << 1) as usize);

                *ptr_y = ia_add_flt(x0r , x1r);
                *(ptr_y.add(1)) = ia_add_flt(x0i, x1i);
                ptr_y = ptr_y.add(2);
            }
        }
    }

    unsafe {
        for i in 0..n_points {
            ptr_real[i as usize] = *y.offset(2 * i as isize);
            ptr_imag[i as usize] = *y.offset((2 * i + 1) as isize);
            //println!("{}  {}", *y.offset(2 * i as isize), *y.offset((2 * i + 1) as isize));
        }
    }  
}


pub const IA_FFT_TWIDDLE_TABLE_FLOAT: [f32; 514] = [
    1.00000000000000000000 ,  0.99998117528260111000 ,  0.99992470183914450000 ,
    0.99983058179582340000 ,  0.99969881869620425000 ,  0.99952941750109314000 ,
    0.99932238458834954000 ,  0.99907772775264536000 ,  0.99879545620517241000 ,
    0.99847558057329477000 ,  0.99811811290014918000 ,  0.99772306664419164000 ,
    0.99729045667869021000 ,  0.99682029929116567000 ,  0.99631261218277800000 ,
    0.99576741446765982000 ,  0.99518472667219693000 ,  0.99456457073425542000 ,
    0.99390697000235606000 ,  0.99321194923479450000 ,  0.99247953459870997000 ,
    0.99170975366909953000 ,  0.99090263542778001000 ,  0.99005821026229712000 ,
    0.98917650996478101000 ,  0.98825756773074946000 ,  0.98730141815785843000 ,
    0.98630809724459867000 ,  0.98527764238894122000 ,  0.98421009238692903000 ,
    0.98310548743121629000 ,  0.98196386910955524000 ,  0.98078528040323043000 ,
    0.97956976568544052000 ,  0.97831737071962765000 ,  0.97702814265775439000 ,
    0.97570213003852857000 ,  0.97433938278557586000 ,  0.97293995220556018000 ,
    0.97150389098625178000 ,  0.97003125319454397000 ,  0.96852209427441738000 ,
    0.96697647104485207000 ,  0.96539444169768940000 ,  0.96377606579543984000 ,
    0.96212140426904158000 ,  0.96043051941556579000 ,  0.95870347489587160000 ,
    0.95694033573220882000 ,  0.95514116830577078000 ,  0.95330604035419386000 ,
    0.95143502096900834000 ,  0.94952818059303667000 ,  0.94758559101774109000 ,
    0.94560732538052128000 ,  0.94359345816196039000 ,  0.94154406518302081000 ,
    0.93945922360218992000 ,  0.93733901191257496000 ,  0.93518350993894761000 ,
    0.93299279883473896000 ,  0.93076696107898371000 ,  0.92850608047321559000 ,
    0.92621024213831138000 ,  0.92387953251128674000 ,  0.92151403934204201000 ,
    0.91911385169005777000 ,  0.91667905992104270000 ,  0.91420975570353069000 ,
    0.91170603200542988000 ,  0.90916798309052238000 ,  0.90659570451491533000 ,
    0.90398929312344334000 ,  0.90134884704602203000 ,  0.89867446569395382000 ,
    0.89596624975618522000 ,  0.89322430119551532000 ,  0.89044872324475788000 ,
    0.88763962040285393000 ,  0.88479709843093779000 ,  0.88192126434835505000 ,
    0.87901222642863353000 ,  0.87607009419540660000 ,  0.87309497841829009000 ,
    0.87008699110871146000 ,  0.86704624551569265000 ,  0.86397285612158681000 ,
    0.86086693863776731000 ,  0.85772861000027212000 ,  0.85455798836540053000 ,
    0.85135519310526520000 ,  0.84812034480329723000 ,  0.84485356524970712000 ,
    0.84155497743689844000 ,  0.83822470555483808000 ,  0.83486287498638001000 ,
    0.83146961230254524000 ,  0.82804504525775580000 ,  0.82458930278502529000 ,
    0.82110251499110465000 ,  0.81758481315158371000 ,  0.81403632970594841000 ,
    0.81045719825259477000 ,  0.80684755354379933000 ,  0.80320753148064494000 ,
    0.79953726910790501000 ,  0.79583690460888357000 ,  0.79210657730021239000 ,
    0.78834642762660634000 ,  0.78455659715557524000 ,  0.78073722857209449000 ,
    0.77688846567323244000 ,  0.77301045336273699000 ,  0.76910333764557970000 ,
    0.76516726562245896000 ,  0.76120238548426178000 ,  0.75720884650648457000 ,
    0.75318679904361252000 ,  0.74913639452345937000 ,  0.74505778544146606000 ,
    0.74095112535495911000 ,  0.73681656887736990000 ,  0.73265427167241282000 ,
    0.72846439044822520000 ,  0.72424708295146700000 ,  0.72000250796138165000 ,
    0.71573082528381859000 ,  0.71143219574521643000 ,  0.70710678118654757000 ,
    0.70275474445722530000 ,  0.69837624940897292000 ,  0.69397146088965400000 ,
    0.68954054473706694000 ,  0.68508366777270036000 ,  0.68060099779545313000 ,
    0.67609270357531603000 ,  0.67155895484701833000 ,  0.66699992230363747000 ,
    0.66241577759017178000 ,  0.65780669329707864000 ,  0.65317284295377676000 ,
    0.64851440102211255000 ,  0.64383154288979150000 ,  0.63912444486377573000 ,
    0.63439328416364549000 ,  0.62963823891492710000 ,  0.62485948814238645000 ,
    0.62005721176328921000 ,  0.61523159058062682000 ,  0.61038280627630948000 ,
    0.60551104140432555000 ,  0.60061647938386897000 ,  0.59569930449243347000 ,
    0.59075970185887428000 ,  0.58579785745643886000 ,  0.58081395809576453000 ,
    0.57580819141784534000 ,  0.57078074588696737000 ,  0.56573181078361323000 ,
    0.56066157619733603000 ,  0.55557023301960229000 ,  0.55045797293660481000 ,
    0.54532498842204646000 ,  0.54017147272989297000 ,  0.53499761988709726000 ,
    0.52980362468629483000 ,  0.52458968267846884000 ,  0.51935599016558953000 ,
    0.51410274419322166000 ,  0.50883014254310699000 ,  0.50353838372571758000 ,
    0.49822766697278187000 ,  0.49289819222978409000 ,  0.48755016014843605000 ,
    0.48218377207912283000 ,  0.47679923006332225000 ,  0.47139673682599781000 ,
    0.46597649576796613000 ,  0.46053871095824001000 ,  0.45508358712634384000 ,
    0.44961132965460660000 ,  0.44412214457042926000 ,  0.43861623853852771000 ,
    0.43309381885315201000 ,  0.42755509343028220000 ,  0.42200027079979979000 ,
    0.41642956009763732000 ,  0.41084317105790391000 ,  0.40524131400498986000 ,
    0.39962419984564679000 ,  0.39399204006104810000 ,  0.38834504669882630000 ,
    0.38268343236508984000 ,  0.37700741021641831000 ,  0.37131719395183760000 ,
    0.36561299780477396000 ,  0.35989503653498828000 ,  0.35416352542049051000 ,
    0.34841868024943451000 ,  0.34266071731199438000 ,  0.33688985339222005000 ,
    0.33110630575987643000 ,  0.32531029216226298000 ,  0.31950203081601575000 ,
    0.31368174039889157000 ,  0.30784964004153498000 ,  0.30200594931922820000 ,
    0.29615088824362396000 ,  0.29028467725446233000 ,  0.28440753721127182000 ,
    0.27851968938505306000 ,  0.27262135544994898000 ,  0.26671275747489842000 ,
    0.26079411791527557000 ,  0.25486565960451463000 ,  0.24892760574572026000 ,
    0.24298017990326398000 ,  0.23702360599436734000 ,  0.23105810828067128000 ,
    0.22508391135979278000 ,  0.21910124015686977000 ,  0.21311031991609136000 ,
    0.20711137619221856000 ,  0.20110463484209196000 ,  0.19509032201612833000 ,
    0.18906866414980628000 ,  0.18303988795514106000 ,  0.17700422041214886000 ,
    0.17096188876030136000 ,  0.16491312048997009000 ,  0.15885814333386139000 ,
    0.15279718525844341000 ,  0.14673047445536175000 ,  0.14065823933284924000 ,
    0.13458070850712622000 ,  0.12849811079379322000 ,  0.12241067519921628000 ,
    0.11631863091190488000 ,  0.11022220729388318000 ,  0.10412163387205473000 ,
    0.09801714032956077000 ,  0.09190895649713269600 ,  0.08579731234443988000 ,
    0.07968243797143012600 ,  0.07356456359966745400 ,  0.06744391956366410600 ,
    0.06132073630220864800 ,  0.05519524434969003100 ,  0.04906767432741812600 ,
    0.04293825693494095900 ,  0.03680722294135899100 ,  0.03067480317663658100 ,
    0.02454122852291226400 ,  0.01840672990580482000 ,  0.01227153828571994400 ,
    0.00613588464915451520 ,  0.00000000000000006123 ,  0.00000000000000000000 ,
    -0.00613588464915447530 , -0.01227153828571992500 , -0.01840672990580482000 ,
    -0.02454122852291228800 , -0.03067480317663662600 , -0.03680722294135883200 ,
    -0.04293825693494082000 , -0.04906767432741801500 , -0.05519524434968993400 ,
    -0.06132073630220857800 , -0.06744391956366405100 , -0.07356456359966742600 ,
    -0.07968243797143012600 , -0.08579731234443989400 , -0.09190895649713272400 ,
    -0.09801714032956060400 , -0.10412163387205459000 , -0.11022220729388306000 ,
    -0.11631863091190475000 , -0.12241067519921620000 , -0.12849811079379317000 ,
    -0.13458070850712617000 , -0.14065823933284921000 , -0.14673047445536175000 ,
    -0.15279718525844344000 , -0.15885814333386145000 , -0.16491312048996989000 ,
    -0.17096188876030122000 , -0.17700422041214875000 , -0.18303988795514095000 ,
    -0.18906866414980619000 , -0.19509032201612825000 , -0.20110463484209190000 ,
    -0.20711137619221856000 , -0.21311031991609136000 , -0.21910124015686980000 ,
    -0.22508391135979283000 , -0.23105810828067111000 , -0.23702360599436720000 ,
    -0.24298017990326387000 , -0.24892760574572015000 , -0.25486565960451457000 ,
    -0.26079411791527551000 , -0.26671275747489837000 , -0.27262135544994898000 ,
    -0.27851968938505306000 , -0.28440753721127188000 , -0.29028467725446233000 ,
    -0.29615088824362379000 , -0.30200594931922808000 , -0.30784964004153487000 ,
    -0.31368174039889152000 , -0.31950203081601569000 , -0.32531029216226293000 ,
    -0.33110630575987643000 , -0.33688985339222005000 , -0.34266071731199438000 ,
    -0.34841868024943456000 , -0.35416352542049034000 , -0.35989503653498811000 ,
    -0.36561299780477385000 , -0.37131719395183754000 , -0.37700741021641826000 ,
    -0.38268343236508978000 , -0.38834504669882625000 , -0.39399204006104810000 ,
    -0.39962419984564679000 , -0.40524131400498986000 , -0.41084317105790391000 ,
    -0.41642956009763715000 , -0.42200027079979968000 , -0.42755509343028208000 ,
    -0.43309381885315196000 , -0.43861623853852766000 , -0.44412214457042920000 ,
    -0.44961132965460654000 , -0.45508358712634384000 , -0.46053871095824001000 ,
    -0.46597649576796618000 , -0.47139673682599764000 , -0.47679923006332209000 ,
    -0.48218377207912272000 , -0.48755016014843600000 , -0.49289819222978404000 ,
    -0.49822766697278187000 , -0.50353838372571758000 , -0.50883014254310699000 ,
    -0.51410274419322166000 , -0.51935599016558964000 , -0.52458968267846895000 ,
    -0.52980362468629461000 , -0.53499761988709715000 , -0.54017147272989285000 ,
    -0.54532498842204646000 , -0.55045797293660481000 , -0.55557023301960218000 ,
    -0.56066157619733603000 , -0.56573181078361312000 , -0.57078074588696726000 ,
    -0.57580819141784534000 , -0.58081395809576453000 , -0.58579785745643886000 ,
    -0.59075970185887416000 , -0.59569930449243336000 , -0.60061647938386897000 ,
    -0.60551104140432555000 , -0.61038280627630948000 , -0.61523159058062682000 ,
    -0.62005721176328910000 , -0.62485948814238634000 , -0.62963823891492698000 ,
    -0.63439328416364549000 , -0.63912444486377573000 , -0.64383154288979139000 ,
    -0.64851440102211244000 , -0.65317284295377676000 , -0.65780669329707864000 ,
    -0.66241577759017178000 , -0.66699992230363747000 , -0.67155895484701833000 ,
    -0.67609270357531592000 , -0.68060099779545302000 , -0.68508366777270036000 ,
    -0.68954054473706683000 , -0.69397146088965400000 , -0.69837624940897292000 ,
    -0.70275474445722530000 , -0.70710678118654746000 , -0.71143219574521643000 ,
    -0.71573082528381859000 , -0.72000250796138165000 , -0.72424708295146689000 ,
    -0.72846439044822520000 , -0.73265427167241282000 , -0.73681656887736979000 ,
    -0.74095112535495911000 , -0.74505778544146595000 , -0.74913639452345926000 ,
    -0.75318679904361241000 , -0.75720884650648446000 , -0.76120238548426178000 ,
    -0.76516726562245896000 , -0.76910333764557959000 , -0.77301045336273699000 ,
    -0.77688846567323244000 , -0.78073722857209438000 , -0.78455659715557524000 ,
    -0.78834642762660623000 , -0.79210657730021239000 , -0.79583690460888346000 ,
    -0.79953726910790501000 , -0.80320753148064483000 , -0.80684755354379922000 ,
    -0.81045719825259477000 , -0.81403632970594830000 , -0.81758481315158371000 ,
    -0.82110251499110465000 , -0.82458930278502529000 , -0.82804504525775580000 ,
    -0.83146961230254524000 , -0.83486287498638001000 , -0.83822470555483797000 ,
    -0.84155497743689833000 , -0.84485356524970701000 , -0.84812034480329712000 ,
    -0.85135519310526520000 , -0.85455798836540053000 , -0.85772861000027212000 ,
    -0.86086693863776731000 , -0.86397285612158670000 , -0.86704624551569265000 ,
    -0.87008699110871135000 , -0.87309497841829009000 , -0.87607009419540660000 ,
    -0.87901222642863341000 , -0.88192126434835494000 , -0.88479709843093779000 ,
    -0.88763962040285393000 , -0.89044872324475788000 , -0.89322430119551532000 ,
    -0.89596624975618511000 , -0.89867446569395382000 , -0.90134884704602203000 ,
    -0.90398929312344334000 , -0.90659570451491533000 , -0.90916798309052227000 ,
    -0.91170603200542988000 , -0.91420975570353069000 , -0.91667905992104270000 ,
    -0.91911385169005777000 , -0.92151403934204190000 , -0.92387953251128674000 ,
    -0.92621024213831127000 , -0.92850608047321548000 , -0.93076696107898371000 ,
    -0.93299279883473885000 , -0.93518350993894750000 , -0.93733901191257496000 ,
    -0.93945922360218992000 , -0.94154406518302081000 , -0.94359345816196039000 ,
    -0.94560732538052128000 , -0.94758559101774109000 , -0.94952818059303667000 ,
    -0.95143502096900834000 , -0.95330604035419375000 , -0.95514116830577067000 ,
    -0.95694033573220894000 , -0.95870347489587160000 , -0.96043051941556579000 ,
    -0.96212140426904158000 , -0.96377606579543984000 , -0.96539444169768940000 ,
    -0.96697647104485207000 , -0.96852209427441727000 , -0.97003125319454397000 ,
    -0.97150389098625178000 , -0.97293995220556007000 , -0.97433938278557586000 ,
    -0.97570213003852857000 , -0.97702814265775439000 , -0.97831737071962765000 ,
    -0.97956976568544052000 , -0.98078528040323043000 , -0.98196386910955524000 ,
    -0.98310548743121629000 , -0.98421009238692903000 , -0.98527764238894122000 ,
    -0.98630809724459867000 , -0.98730141815785843000 , -0.98825756773074946000 ,
    -0.98917650996478101000 , -0.99005821026229712000 , -0.99090263542778001000 ,
    -0.99170975366909953000 , -0.99247953459870997000 , -0.99321194923479450000 ,
    -0.99390697000235606000 , -0.99456457073425542000 , -0.99518472667219682000 ,
    -0.99576741446765982000 , -0.99631261218277800000 , -0.99682029929116567000 ,
    -0.99729045667869021000 , -0.99772306664419164000 , -0.99811811290014918000 ,
    -0.99847558057329477000 , -0.99879545620517241000 , -0.99907772775264536000 ,
    -0.99932238458834954000 , -0.99952941750109314000 , -0.99969881869620425000 ,
    -0.99983058179582340000 , -0.99992470183914450000 , -0.99998117528260111000 ,
    -1.00000000000000000000];
