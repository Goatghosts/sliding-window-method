// Source: https://github.com/ingonyama-zk/papers/blob/main/modular_multiplication.pdf
`timescale 1ns/1ps
module ModMul #(
  parameter p = 37, // a*b mod p = r
  parameter width = 128 // width of the output
) (
  input clk,
  input reset,
  input [width-1:0] a,        // First multiplication element, a
  input [width-1:0] b,        // Second multiplication element, b
  input enable,
  output [width-1:0] r,         // Remainder, r
  output done
);

  wire mul_done;
  wire [2*width-1:0] ab;
  reg reset_reducer, enable_reducer;

  // Booth only
  reg ld;

  always @(posedge clk) begin
    if (reset) begin
      ld <= 1'b0;
      reset_reducer <= 0;
      enable_reducer <= 0;
    end else begin
      if (ld) begin
        ld <= 1'b0;
      end
      if (mul_done && !enable_reducer) begin
        reset_reducer <= 1;
        enable_reducer <= 1;
      end else begin
        reset_reducer <= 0;
      end
    end
  end
  

  ModReduction #(
    .p(p),
    .width(width)
  )
  u_mod_reduction (
      .clk(clk),
      .reset(reset_reducer),
      .enable(enable_reducer),
      .a(ab),
      .done(done),
      .r(r)
  );

  // Multiplication module
  // - Karatsuba
 karat_mult_recursion #(
   .wI         (width),
   .nSTAGE     (5)
 )
 u_karat_mult_recursion (
   // data IOs
   .iX(a),
   .iY(b),
   .oO(ab),
   // control IOs
   .clk(clk),
   .reset(reset),
   .i_enable(enable),
   .o_finish(mul_done) 
 );

  // - Booth multiplication

// always @(posedge enable) begin
//   ld <= 1'b1;
// end

// Booth_Multiplier #(
//     .pN(7)
//   ) u_booth_mult (
//     .Rst(reset), 
//     .Clk(clk), 
//     .Ld(ld), 
//     .M(a), 
//     .R(b), 
//     .Valid(mul_done), 
//     .P(ab)
//   );

endmodule